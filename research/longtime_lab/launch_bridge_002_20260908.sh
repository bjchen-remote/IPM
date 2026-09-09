#!/bin/zsh
set -eu
campaign_project='/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured'
campaign_root="${campaign_project}/result/longtime/20260908_campaign_v1"
test ! -d "${campaign_root}/adaptive.lock"
mkdir "${campaign_root}/bridge.lock"
trap 'rmdir "${campaign_root}/bridge.lock" 2>/dev/null || true' EXIT
cd "${campaign_project}"
/Applications/MATLAB_R2026a.app/bin/matlab -batch "run('result/longtime/20260908_campaign_v1/bridge_run_002.m')" > "${campaign_root}/bridge_002.log" 2>&1
