#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
release_parent="${1:-$(dirname "$project_root")/releases}"
release_name="${2:-ipm_long_time_server_v5_open_horizon_20260913}"
if [[ ! "$release_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "Release name must be a single safe directory name." >&2
    exit 2
fi
stage="$release_parent/$release_name"
archive="$release_parent/$release_name.zip"

mkdir -p "$release_parent"
if [[ -e "$stage" || -e "$archive" || -e "$archive.sha256" ]]; then
    echo "Refusing to overwrite an existing release artifact." >&2
    exit 2
fi
mkdir "$stage"
cp -R "$project_root/+ipm" "$stage/+ipm"
cp -R "$project_root/tests" "$stage/tests"
cp -R "$project_root/examples" "$stage/examples"
cp -R "$project_root/server" "$stage/server"
cp "$project_root/README.md" "$stage/README.md"
cp "$project_root/README_SERVER_ZH.md" "$stage/README_SERVER_ZH.md"
cp "$project_root/MESH_GRID_RECOMMENDATION_ZH.md" "$stage/MESH_GRID_RECOMMENDATION_ZH.md"
cp "$project_root/STRUCTURE.md" "$stage/STRUCTURE.md"
cp "$project_root/CHANGELOG.md" "$stage/CHANGELOG.md"
mkdir -p "$stage/research/longtime_lab/evidence"
cp "$project_root/research/longtime_lab/AUTONOMOUS_GRID_DESIGN_FROM_FROZEN_DATA_20260913.md" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/AMORTIZED_MESH_LIFETIME_HOLDOUT_20260913.md" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/CROSS_LEVEL_MESH_COST_AND_RHS_20260913.md" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/evidence/"*.json "$stage/research/longtime_lab/evidence/"
cp "$project_root/research/longtime_lab/audit_frozen_pair_ratio.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/verify_continuous_vertical_shadow.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/run_shadow_resume_smoke.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/verify_shadow_native_equivalence.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/audit_amortized_mesh_ratio.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/audit_alternative_candidate_transfer.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/replay_requested_alternative_transfer.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/rollout_alternative_mesh_interval.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/audit_candidate_lifetime_at_request.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/audit_cross_level_mesh_at_request.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/rollout_cross_level_mesh_interval.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/verify_ranked_lifetime_followup.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/continuous_inner_geometry.m" "$stage/research/longtime_lab/"
cp "$project_root/research/longtime_lab/render_continuous_profile_progress.m" "$stage/research/longtime_lab/"
mkdir -p "$stage/research/acceleration_lab/evidence"
cp "$project_root/research/acceleration_lab/ORIGINAL_T0_TAU13_SHAPE_AND_REMESH_20260913.md" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ORIGINAL_T0_TAU139_PROFILE_20260913.md" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/REMESH_RHS_C1_DECOMPOSITION_20260913.md" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/audit_native_remesh_profile_jump.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/audit_native_remesh_shape_tangent.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/audit_original_t0_shape_train.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/audit_native_remesh_rhs_c_decomposition.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ipm_accellab_continuous_inner_rates.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ipm_accellab_hermite_level_roots.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ipm_accellab_hermite_line.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ipm_accellab_hermite_peak.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/ipm_accellab_tensor_hermite.m" "$stage/research/acceleration_lab/"
cp "$project_root/research/acceleration_lab/evidence/"*.json "$stage/research/acceleration_lab/evidence/"
cp "$project_root/research/acceleration_lab/evidence/"*.png "$stage/research/acceleration_lab/evidence/"
cp "$project_root/research/acceleration_lab/evidence/"*.mat "$stage/research/acceleration_lab/evidence/"
chmod +x "$stage/server/launch.sh" "$stage/server/package_release.sh"

commit="$(git -C "$project_root" rev-parse HEAD)"
if [[ -z "$(git -C "$project_root" status --porcelain)" ]]; then
    working_tree_clean=true
else
    working_tree_clean=false
fi
{
    echo "release=$release_name"
    echo "created_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "source_commit=$commit"
    echo "working_tree_clean_at_packaging=$working_tree_clean"
    echo "matlab_reference=R2026a"
} > "$stage/PACKAGE_INFO.txt"
(
    cd "$stage"
    find . -type f ! -name SHA256SUMS -print | LC_ALL=C sort | \
        while IFS= read -r file; do shasum -a 256 "$file"; done > SHA256SUMS
)
(
    cd "$release_parent"
    zip -qr "$archive" "$release_name"
)
shasum -a 256 "$archive" > "$archive.sha256"
echo "$archive"
