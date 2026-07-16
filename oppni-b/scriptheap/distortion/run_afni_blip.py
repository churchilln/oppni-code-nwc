#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from pathlib import Path


def run(cmd):
    """Run one shell command and stop immediately if it fails."""
    print("\n>>>", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


def require_tool(name):
    """Make sure the needed AFNI command is available before starting."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Cannot find command: {name}")


def main():
    parser = argparse.ArgumentParser(
        description="Run AFNI AP/PA blip distortion correction."
    )

    parser.add_argument("--ap", required=True, help="Main AP BOLD NIfTI")
    parser.add_argument("--pa", required=True, help="Reverse PA BOLD NIfTI")
    parser.add_argument("--subj-id", default="sub039_afni_blip", help="AFNI subject ID")
    parser.add_argument("--outdir", default="afni_blip_output", help="Output folder")

    args = parser.parse_args()

    # Fail early if AFNI is not properly loaded.
    for tool in ["afni_proc.py", "3dQwarp", "tcsh", "3dTstat"]:
        require_tool(tool)

    ap = Path(args.ap).resolve()
    pa = Path(args.pa).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    proc_script = outdir / f"proc.{args.subj_id}"
    results_dir = outdir / f"{args.subj_id}.results"

    # Run AFNI AP/PA blip correction with basic fMRI preprocessing blocks.
    run([
        "afni_proc.py",
        "-subj_id", args.subj_id,
        "-script", proc_script,
        "-out_dir", results_dir,
        "-dsets", ap,
        "-blip_forward_dset", ap,
        "-blip_reverse_dset", pa,
        "-blip_opts_qw", "-noXdis", "-noZdis",
        "-blocks", "tshift", "volreg",
        "-volreg_align_to", "MIN_OUTLIER",
        "-html_review_style", "pythonic",
        "-execute",
    ])

    afni_blip_output = results_dir / f"pb02.{args.subj_id}.r01.blip+orig"
    afni_blip_mean = results_dir / "AFNI_blip_corrected_mean.nii.gz"

    # Mean AFNI blip-corrected AP is only for visual comparison.
    run([
        "3dTstat",
        "-mean",
        "-prefix", afni_blip_mean,
        afni_blip_output,
    ])

    print("\nDONE")
    print(f"AFNI script: {proc_script}")
    print(f"Results folder: {results_dir}")
    print(f"AFNI blip-corrected AP: {afni_blip_output}")
    print(f"Corrected mean: {afni_blip_mean}")
    print("\nCheck the corrected mean and QC output.")


if __name__ == "__main__":
    main()
