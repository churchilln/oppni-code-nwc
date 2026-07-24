#!/usr/bin/env python3
"""Warp OPPNI-D scalar maps with TBSS-cleaned FA and ANTs."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


DTI_SCALAR_RE = re.compile(r"^dtifit_[0-9]+_(FA|MD|L1|L2|L3|S0)\.nii(\.gz)?$")
DKI_SCALAR_NAMES = (
    "dtifit_ms_kurt.nii.gz",
    "dtifit_ms_kurt1.nii.gz",
    "dtifit_ms_kurt2.nii.gz",
    "dtifit_ms_kurt3.nii.gz",
    "dtifit_ms_S0.nii.gz",
)


# parse subject and template inputs supplied by MATLAB
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Warp OPPNI-D DTI/DKI scalar maps.")
    parser.add_argument("--subject-dir", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


# locate required external commands
def find_command(command: str) -> str:
    candidates = [shutil.which(command)]
    candidates.append(str(Path(sys.executable).parent / command))
    fsldir = os.environ.get("FSLDIR")
    if fsldir:
        candidates.append(str(Path(fsldir) / "bin" / command))

    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return candidate

    raise RuntimeError(f"cannot find required command: {command}")


# run command and save stdout/stderr to log file
def run_command(command: list[str], log_file: Path, cwd: Path | None, env: dict[str, str]) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    printable_command = " ".join(command)
    print(f"\nrunning:\n{printable_command}\n")

    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    with log_file.open("a", encoding="utf-8") as log:
        log.write("\n")
        log.write("=" * 80)
        log.write(f"\nCOMMAND:\n{printable_command}\n\n")
        log.write(result.stdout)
        log.write("\n")

    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"command failed with exit code {result.returncode}: {printable_command}")


def nii_stem(image: Path) -> str:
    if image.name.endswith(".nii.gz"):
        return image.name[:-7]
    if image.name.endswith(".nii"):
        return image.name[:-4]
    return image.stem


# select DTI scalar maps; exclude vectors/tensors
def find_dti_scalars(dti_dir: Path) -> list[Path]:
    if not dti_dir.is_dir():
        return []
    return [path.resolve() for path in sorted(dti_dir.glob("*.nii*")) if DTI_SCALAR_RE.match(path.name)]


# select DKI scalar maps; exclude tensor output
def find_dki_scalars(dki_dir: Path) -> list[Path]:
    if not dki_dir.is_dir():
        return []
    return [path.resolve() for name in DKI_SCALAR_NAMES for path in sorted(dki_dir.glob(name))]


def make_clean_fa(
    native_fa: Path,
    subject_id: str,
    temp_tbss_dir: Path,
    log_file: Path,
    env: dict[str, str],
    tbss_preproc: str,
) -> Path:
    # run TBSS preprocessing to generate registration FA
    if temp_tbss_dir.exists():
        shutil.rmtree(temp_tbss_dir)
    temp_tbss_dir.mkdir(parents=True)

    tbss_input = temp_tbss_dir / f"{subject_id}_FA.nii.gz"
    shutil.copy2(native_fa, tbss_input)

    run_command([tbss_preproc, tbss_input.name], log_file, temp_tbss_dir, env)

    clean_fa = temp_tbss_dir / "FA" / f"{subject_id}_FA_FA.nii.gz"
    if not clean_fa.is_file():
        raise RuntimeError(f"expected TBSS-cleaned FA was not created: {clean_fa}")
    return clean_fa


def apply_ants_transform(
    input_image: Path,
    output_image: Path,
    template: Path,
    warp: Path,
    affine: Path,
    log_file: Path,
    env: dict[str, str],
    ants_apply: str,
) -> None:
    output_image.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            ants_apply,
            "-d",
            "3",
            "-i",
            str(input_image),
            "-r",
            str(template),
            "-t",
            str(warp),
            "-t",
            str(affine),
            "-n",
            "Linear",
            "-o",
            str(output_image),
        ],
        log_file,
        None,
        env,
    )
    if not output_image.is_file():
        raise RuntimeError(f"expected warped image was not created: {output_image}")


def main() -> int:
    args = parse_args()
    subject_dir = args.subject_dir.expanduser().resolve()
    template = args.template.expanduser().resolve()
    subject_id = subject_dir.name

    tbss_preproc = find_command("tbss_1_preproc")
    ants_registration = find_command("antsRegistrationSyN.sh")
    ants_apply = find_command("antsApplyTransforms")

    p2_dir = subject_dir / "diff_proc_p2"
    dti_dir = p2_dir / "dti"
    dki_dir = p2_dir / "dki"
    native_fa = dti_dir / "dtifit_1_FA.nii.gz"

    if not p2_dir.is_dir():
        raise FileNotFoundError(f"diff_proc_p2 directory not found: {p2_dir}")
    if not native_fa.is_file():
        raise FileNotFoundError(f"registration FA not found: {native_fa}")
    if not template.is_file():
        raise FileNotFoundError(f"FA template not found: {template}")

    align_dir = p2_dir / "alignment" / "ants_tbss"
    temp_tbss_dir = align_dir / "temp_tbss"
    transform_dir = align_dir / "transforms"
    qc_dir = align_dir / "qc"
    log_file = align_dir / "logs" / "registration.log"
    status_file = align_dir / "status.json"

    transform_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(max(1, args.threads))

    clean_fa = make_clean_fa(native_fa, subject_id, temp_tbss_dir, log_file, env, tbss_preproc)
    shutil.copy2(clean_fa, qc_dir / "registration_FA_cleaned_native.nii.gz")

    ants_prefix = transform_dir / "FA_to_FMRIB58_"
    run_command(
        [
            ants_registration,
            "-d",
            "3",
            "-f",
            str(template),
            "-m",
            str(clean_fa),
            "-o",
            str(ants_prefix),
            "-t",
            "s",
        ],
        log_file,
        None,
        env,
    )

    affine = Path(f"{ants_prefix}0GenericAffine.mat")
    warp = Path(f"{ants_prefix}1Warp.nii.gz")
    inverse_warp = Path(f"{ants_prefix}1InverseWarp.nii.gz")
    warped_fa = Path(f"{ants_prefix}Warped.nii.gz")
    for path in (affine, warp, inverse_warp, warped_fa):
        if not path.is_file():
            raise RuntimeError(f"expected ANTs output was not created: {path}")
    shutil.copy2(warped_fa, qc_dir / "registration_FA_cleaned_FMRIB58.nii.gz")

    # apply the FA-derived transform to scalar model outputs
    template_label = nii_stem(template)
    warped_outputs: dict[str, list[str]] = {"dti": [], "dki": []}

    for input_image in find_dti_scalars(dti_dir):
        output_image = dti_dir / "warped" / f"{nii_stem(input_image)}_{template_label}.nii.gz"
        apply_ants_transform(input_image, output_image, template, warp, affine, log_file, env, ants_apply)
        warped_outputs["dti"].append(str(output_image))

    for input_image in find_dki_scalars(dki_dir):
        output_image = dki_dir / "warped" / f"{nii_stem(input_image)}_{template_label}.nii.gz"
        apply_ants_transform(input_image, output_image, template, warp, affine, log_file, env, ants_apply)
        warped_outputs["dki"].append(str(output_image))

    if not warped_outputs["dti"]:
        raise RuntimeError(f"no DTI scalar maps were found in: {dti_dir}")

    status = {
        "status": "complete",
        "subject_dir": str(subject_dir),
        "registration_fa": str(native_fa),
        "tbss_cleaned_fa": str(clean_fa),
        "template": str(template),
        "transforms": {
            "affine": str(affine),
            "warp": str(warp),
            "inverse_warp": str(inverse_warp),
        },
        "warped_outputs": warped_outputs,
    }
    with status_file.open("w", encoding="utf-8") as file_obj:
        json.dump(status, file_obj, indent=2)

    print("\ndiffusion warp complete")
    print(f"alignment: {align_dir}")
    print(f"DTI warped outputs: {dti_dir / 'warped'}")
    print(f"DKI warped outputs: {dki_dir / 'warped'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
