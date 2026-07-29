#!/usr/bin/env bash
set -euo pipefail

OPPNI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$OPPNI_ROOT/.venv"
REQ="$OPPNI_ROOT/requirements.txt"
PYTHON_BIN="${PYTHON:-python3}"

unset PYTHONPATH

if [[ ! -x "$VENV/bin/python" ]]; then
    "$PYTHON_BIN" -m venv "$VENV"
fi

# Install a CUDA-compatible PyTorch build on Compute Canada clusters.
if [[ -d /cvmfs/soft.computecanada.ca/custom/python/wheelhouse ]]; then
    "$VENV/bin/python" -m pip install "torch==2.6.0+computecanada"
else
    "$VENV/bin/python" -m pip install "torch==2.6.0"
fi

"$VENV/bin/python" -m pip install -r "$REQ"
"$VENV/bin/python" -c "from deepbet import run_bet"

echo "OPPNI Python setup complete."
