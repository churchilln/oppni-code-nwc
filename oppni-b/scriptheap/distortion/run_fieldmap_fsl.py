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
    """Load JSON metadata."""
    with open(path, "r") as f:
        return json.load(f)


def require_tool(name):
    """Make sure the needed FSL command is available before starting."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Cannot find FSL command: {name}")


def pe_to_unwarpdir(pe):
    """Convert BIDS phase-encoding direction to FSL FUGUE unwarpdir format."""
    mapping = {
        "i": "x",
        "i-": "x-",
        "j": "y",
        "j-": "y-",
        "k": "z",
        "k-": "z-",
    }

    if pe not in mapping:
        raise ValueError(f"Unknown PhaseEncodingDirection: {pe}")

    return mapping[pe]


def require_value(value, label):
    """Return a required command-line or JSON value."""
    if value is None or value == "":
        raise ValueError(f"Missing required value: {label}")
    return value


def get_delta_te_ms(args, fmap_meta):
    """Get fieldmap echo-time difference in milliseconds."""
    if args.echo_time_difference is not None:
        return abs(float(args.echo_time_difference)) * 1000

    if args.echo_time_1 is not None and args.echo_time_2 is not None:
        return abs(float(args.echo_time_2) - float(args.echo_time_1)) * 1000

    if "EchoTime1" in fmap_meta and "EchoTime2" in fmap_meta:
        return abs(fmap_meta["EchoTime2"] - fmap_meta["EchoTime1"]) * 1000

    if "EchoTimeDifference" in fmap_meta:
        return fmap_meta["EchoTimeDifference"] * 1000

    raise ValueError("No EchoTime1/EchoTime2 or EchoTimeDifference found")


def main():
    parser = argparse.ArgumentParser(
        description="Run FSL fieldmap-based distortion correction."
    )

    parser.add_argument("--ap", required=True, help="Main AP BOLD NIfTI")
    parser.add_argument("--ap-json", help="Optional AP BOLD JSON fallback")
    parser.add_argument("--mag", required=True, help="Fieldmap magnitude NIfTI, usually magnitude1")
    parser.add_argument("--phasediff", required=True, help="Fieldmap phasediff NIfTI")
    parser.add_argument("--phasediff-json", help="Optional fieldmap phasediff JSON fallback")
    parser.add_argument("--phase-encoding-direction", help="AP BOLD PhaseEncodingDirection")
    parser.add_argument("--effective-echo-spacing", help="AP BOLD EffectiveEchoSpacing")
    parser.add_argument("--echo-time-1", help="Fieldmap first echo time in seconds")
    parser.add_argument("--echo-time-2", help="Fieldmap second echo time in seconds")
    parser.add_argument("--echo-time-difference", help="Fieldmap echo-time difference in seconds")
    parser.add_argument("--outdir", default="fieldmap_fsl_output", help="Output folder")

    args = parser.parse_args()

    # Fail early if FSL is not properly loaded.
    for tool in ["fslmaths", "bet", "fsl_prepare_fieldmap", "fugue"]:
        require_tool(tool)

    ap = Path(args.ap).resolve()
    mag = Path(args.mag).resolve()
    phasediff = Path(args.phasediff).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # MATLAB passes explicit values; JSON remains as a standalone fallback.
    ap_meta = read_json(Path(args.ap_json).resolve()) if args.ap_json else {}
    fmap_meta = read_json(Path(args.phasediff_json).resolve()) if args.phasediff_json else {}

    ap_pe = require_value(args.phase_encoding_direction or ap_meta.get("PhaseEncodingDirection"), "AP phase-encoding direction")
    unwarpdir = pe_to_unwarpdir(ap_pe)
    dwell = require_value(args.effective_echo_spacing or ap_meta.get("EffectiveEchoSpacing"), "effective echo spacing")
    delta_te_ms = get_delta_te_ms(args, fmap_meta)

    print(f"AP PhaseEncodingDirection: {ap_pe}")
    print(f"FSL unwarpdir: {unwarpdir}")
    print(f"AP EffectiveEchoSpacing: {dwell}")
    print(f"Fieldmap Delta TE: {delta_te_ms:.6f} ms")

    fmap_mag_brain = outdir / "fmap_mag_brain.nii.gz"
    fmap_rads = outdir / "fmap_rads.nii.gz"

    ap_uncorrected_mean = outdir / "AP_uncorrected_mean.nii.gz"
    ap_corrected = outdir / "AP_fieldmap_corrected.nii.gz"
    ap_corrected_mean = outdir / "AP_fieldmap_corrected_mean.nii.gz"

    # Mean original AP is only for visual comparison.
    run(["fslmaths", ap, "-Tmean", ap_uncorrected_mean])

    # Brain extract the magnitude image for fieldmap preparation.
    run(["bet", mag, fmap_mag_brain, "-f", "0.35", "-R"])

    # Convert Siemens phasediff into an FSL-ready fieldmap in rad/s.
    run([
        "fsl_prepare_fieldmap",
        "SIEMENS",
        phasediff,
        fmap_mag_brain,
        fmap_rads,
        f"{delta_te_ms:.6f}",
    ])

    # Apply the prepared fieldmap to the full AP BOLD.
    run([
        "fugue",
        "-i", ap,
        f"--loadfmap={fmap_rads}",
        f"--dwell={dwell}",
        f"--unwarpdir={unwarpdir}",
        "-u", ap_corrected,
    ])

    # Mean corrected AP is only for visual comparison.
    run(["fslmaths", ap_corrected, "-Tmean", ap_corrected_mean])

    print("\nDONE")
    print(f"Prepared fieldmap: {fmap_rads}")
    print(f"Corrected AP: {ap_corrected}")
    print(f"Uncorrected mean: {ap_uncorrected_mean}")
    print(f"Corrected mean: {ap_corrected_mean}")


if __name__ == "__main__":
    main()
