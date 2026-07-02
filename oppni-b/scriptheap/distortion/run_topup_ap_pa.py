#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path


def run(cmd):
    """Run one shell command and stop immediately if it fails."""
    print("\n>>>", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


def read_json(path):
    """Load BIDS sidecar metadata."""
    with open(path, "r") as f:
        return json.load(f)


def pe_to_vector(pe):
    """Convert BIDS phase-encoding direction to FSL TOPUP vector format."""
    mapping = {
        "i": "1 0 0",
        "i-": "-1 0 0",
        "j": "0 1 0",
        "j-": "0 -1 0",
        "k": "0 0 1",
        "k-": "0 0 -1",
    }

    if pe not in mapping:
        raise ValueError(f"Unknown PhaseEncodingDirection: {pe}")

    return mapping[pe]


def require_tool(name):
    """Make sure the needed FSL command is available before starting."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Cannot find FSL command: {name}")


def require_value(value, label):
    """Return a required command-line or JSON value."""
    if value is None or value == "":
        raise ValueError(f"Missing required value: {label}")
    return value


def main():
    parser = argparse.ArgumentParser(
        description="Run FSL TOPUP AP/PA distortion correction."
    )

    parser.add_argument("--ap", required=True, help="Main AP BOLD NIfTI")
    parser.add_argument("--pa", required=True, help="Reverse PA BOLD NIfTI")
    parser.add_argument("--ap-json", help="Optional AP BOLD JSON fallback")
    parser.add_argument("--pa-json", help="Optional PA BOLD JSON fallback")
    parser.add_argument("--ap-pe-dir", help="AP BOLD PhaseEncodingDirection")
    parser.add_argument("--pa-pe-dir", help="PA BOLD PhaseEncodingDirection")
    parser.add_argument("--ap-total-readout-time", help="AP BOLD TotalReadoutTime")
    parser.add_argument("--pa-total-readout-time", help="PA BOLD TotalReadoutTime")
    parser.add_argument("--outdir", default="topup_ap_pa_output", help="Output folder")

    args = parser.parse_args()

    # Fail early if FSL is not properly loaded.
    for tool in ["fslroi", "fslmerge", "topup", "applytopup", "fslmaths"]:
        require_tool(tool)

    ap = Path(args.ap).resolve()
    pa = Path(args.pa).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # MATLAB passes explicit values; JSON remains as a standalone fallback.
    ap_meta = read_json(Path(args.ap_json).resolve()) if args.ap_json else {}
    pa_meta = read_json(Path(args.pa_json).resolve()) if args.pa_json else {}

    ap_pe = require_value(args.ap_pe_dir or ap_meta.get("PhaseEncodingDirection"), "AP phase-encoding direction")
    pa_pe = require_value(args.pa_pe_dir or pa_meta.get("PhaseEncodingDirection"), "PA phase-encoding direction")
    ap_readout = require_value(args.ap_total_readout_time or ap_meta.get("TotalReadoutTime"), "AP total readout time")
    pa_readout = require_value(args.pa_total_readout_time or pa_meta.get("TotalReadoutTime"), "PA total readout time")

    print(f"AP PhaseEncodingDirection: {ap_pe}")
    print(f"PA PhaseEncodingDirection: {pa_pe}")
    print(f"AP TotalReadoutTime: {ap_readout}")
    print(f"PA TotalReadoutTime: {pa_readout}")

    # TOPUP needs one row per input image, in the same order as fslmerge below.
    acqparams = outdir / "acqparams.txt"

    with open(acqparams, "w") as f:
        f.write(f"{pe_to_vector(ap_pe)} {ap_readout}\n")
        f.write(f"{pe_to_vector(pa_pe)} {pa_readout}\n")

    print(f"\nWrote: {acqparams}")
    print(acqparams.read_text())

    ap_ref = outdir / "AP_ref.nii.gz"
    pa_ref = outdir / "PA_ref.nii.gz"
    blip_pair = outdir / "blip_pair.nii.gz"

    # Extract one volume from AP and PA to estimate the distortion field.
    run(["fslroi", ap, ap_ref, "0", "1"])
    run(["fslroi", pa, pa_ref, "0", "1"])

    # Merge AP first and PA second, matching acqparams row 1 and row 2.
    run(["fslmerge", "-t", blip_pair, ap_ref, pa_ref])

    topup_base = outdir / "topup_out"
    topup_field = outdir / "topup_field.nii.gz"
    topup_corrected_pair = outdir / "topup_corrected_pair.nii.gz"

    # Estimate the distortion field from the AP/PA pair.
    run([
        "topup",
        f"--imain={blip_pair}",
        f"--datain={acqparams}",
        "--config=b02b0.cnf",
        f"--out={topup_base}",
        f"--fout={topup_field}",
        f"--iout={topup_corrected_pair}",
    ])

    ap_corrected = outdir / "AP_topup_corrected.nii.gz"

    # Apply the estimated field to the full AP BOLD.
    # inindex=1 because AP is row 1 in acqparams.txt.
    run([
        "applytopup",
        f"--imain={ap}",
        f"--datain={acqparams}",
        "--inindex=1",
        f"--topup={topup_base}",
        "--method=jac",
        f"--out={ap_corrected}",
    ])

    ap_uncorrected_mean = outdir / "AP_uncorrected_mean.nii.gz"
    ap_corrected_mean = outdir / "AP_topup_corrected_mean.nii.gz"

    # Mean images make visual inspection easier.
    run(["fslmaths", ap, "-Tmean", ap_uncorrected_mean])
    run(["fslmaths", ap_corrected, "-Tmean", ap_corrected_mean])

    print("\nDONE")
    print(f"Corrected AP: {ap_corrected}")
    print(f"Uncorrected mean: {ap_uncorrected_mean}")
    print(f"Corrected mean: {ap_corrected_mean}")


if __name__ == "__main__":
    main()
