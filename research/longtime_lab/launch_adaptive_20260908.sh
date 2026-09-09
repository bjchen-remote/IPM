#!/bin/zsh
set -eu
campaign_project='/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured'
campaign_root="${campaign_project}/result/longtime/20260908_campaign_v1"
campaign_lock="${campaign_root}/adaptive.lock"
mkdir "${campaign_lock}" || exit 2
trap 'rmdir "${campaign_lock}" 2>/dev/null || true' EXIT
export MATLAB_PREFDIR='/private/tmp/ipm-longtime-adaptive-20260908-prefs'
cd "${campaign_project}"
/Applications/MATLAB_R2026a.app/bin/matlab -batch "p=pwd; d=fullfile(p,'result','longtime','20260908_campaign_v1'); s=fullfile(d,'adaptive_source'); addpath(s); addpath(fullfile(s,'research','grid_lab')); addpath(fullfile(s,'research','longtime_lab')); addpath(fullfile(s,'research','gauge_lab')); addpath(fullfile(s,'research','experiments')); cd(s); assert(startsWith(which('ipm.solve'),s)); run_adaptive_campaign(p,fullfile(d,'mesh_tau407_smoke_v1','mesh_smoke_report.mat'),fullfile(d,'adaptive_campaign'),80);" >> "${campaign_root}/adaptive.log" 2>&1
