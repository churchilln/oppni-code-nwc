#!/usr/bin/env bash
set -euo pipefail

OPPNI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$OPPNI_ROOT/.venv"
REQ="$OPPNI_ROOT/requirements.txt"
PYTHON_BIN="${PYTHON:-python3}"

if [[ ! -x "$VENV/bin/python" ]]; then
    "$PYTHON_BIN" -m venv "$VENV"
fi

"$VENV/bin/python" -m pip install -r "$REQ"
"$VENV/bin/python" -c "from deepbet import run_bet"
