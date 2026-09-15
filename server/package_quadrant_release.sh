#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
release_parent="${1:-$(dirname "$project_root")/releases}"
release_name="ipm_long_time_server_20260915"
if [[ $# -gt 1 ]]; then
    echo "Only the release parent may be overridden; main.m fixes the server release name." >&2
    exit 2
fi
expected_main_line="releaseRoot='/data/user/hd58131/ipm/$release_name';"
if ! grep -Fxq "$expected_main_line" "$project_root/server/main_quadrant_release.m"; then
    echo "main.m absolute releaseRoot does not match the package version." >&2
    exit 2
fi
stage="$release_parent/$release_name"
archive="$release_parent/$release_name.zip"
if [[ -e "$stage" || -e "$archive" || -e "$archive.sha256" ]]; then
    echo "Refusing to overwrite an existing quadrant release." >&2
    exit 2
fi

mkdir -p "$release_parent"
mkdir "$stage"
cp -R "$project_root/+ipm" "$stage/+ipm"
cp -R "$project_root/tests" "$stage/tests"
cp -R "$project_root/examples" "$stage/examples"
cp -R "$project_root/server" "$stage/server"
cp -R "$project_root/research" "$stage/research"
cp "$project_root/server/main_quadrant_release.m" "$stage/main.m"
cp "$project_root"/README*.md "$stage/"
cp "$project_root"/RELEASE_NOTES*.md "$stage/"
cp "$project_root/MESH_GRID_RECOMMENDATION_ZH.md" "$stage/"
cp "$project_root/STRUCTURE.md" "$stage/"
cp "$project_root/ARCHITECTURE_QUADRANT_LEVELSET_20260914.md" "$stage/"
cp "$project_root/CHANGELOG.md" "$stage/"

commit="$(git -C "$project_root" rev-parse HEAD)"
if [[ -z "$(git -C "$project_root" status --porcelain)" ]]; then
    working_tree_clean=true
else
    working_tree_clean=false
fi
{
    echo "release=$release_name"
    echo "launcher=quadrant_positive_only_copied_main_v1"
    echo "server_release_root=/data/user/hd58131/ipm/$release_name"
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
(
    cd "$release_parent"
    shasum -a 256 "$release_name.zip" > "$release_name.zip.sha256"
)
echo "$archive"
