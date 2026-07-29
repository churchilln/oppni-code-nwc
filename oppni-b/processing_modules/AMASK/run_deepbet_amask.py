#!/usr/bin/env python3
"""Run DeepBET for OPPNI AMASK."""

import argparse

from deepbet import run_bet


def main():
    parser = argparse.ArgumentParser(description="Create a DeepBET brain and mask.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--brain", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--n-dilate", type=int, default=0)
    args = parser.parse_args()

    run_bet(
        input_paths=[args.input],
        brain_paths=[args.brain],
        mask_paths=[args.mask],
        tiv_paths=None,
        threshold=args.threshold,
        n_dilate=args.n_dilate,
        no_gpu=True,
    )


if __name__ == "__main__":
    main()
