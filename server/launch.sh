#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
matlab_bin="${MATLAB_BIN:-matlab}"
cd "$project_root"
exec "$matlab_bin" -batch "run('server/run_profile.m')"
