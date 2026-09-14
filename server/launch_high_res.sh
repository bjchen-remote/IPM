#!/usr/bin/env bash
set -euo pipefail

release_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export IPM_MAIN_PATH="$release_root/main.m"
matlab_bin="${MATLAB_BIN:-matlab}"
exec "$matlab_bin" -batch "run(getenv('IPM_MAIN_PATH'))"
