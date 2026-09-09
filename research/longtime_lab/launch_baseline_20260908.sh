#!/bin/zsh
set -eu
campaign_project='/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured'
campaign_directory="${campaign_project}/result/longtime/20260908_campaign_v1"
campaign_lock="${campaign_directory}/baseline.lock"
mkdir "${campaign_lock}" || exit 2
trap 'rmdir "${campaign_lock}" 2>/dev/null || true' EXIT
export MATLAB_PREFDIR='/private/tmp/ipm-longtime-baseline-20260908-prefs'
cd "${campaign_project}"
/Applications/MATLAB_R2026a.app/bin/matlab -batch "run('research/longtime_lab/run_baseline_20260908.m')" >> "${campaign_directory}/baseline.log" 2>&1
