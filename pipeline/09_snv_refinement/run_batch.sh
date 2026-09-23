#!/bin/bash
# Run refine_pair.sh on several pairs in series; one failure does not stop the rest.
#
#   bash scripts/run_batch.sh 202 216 226 236 240 246 254 275 279
set -uo pipefail
[ "$#" -gt 0 ] || { echo "usage: bash run_batch.sh <pair_id>..."; exit 1; }
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
failed=()
for id in "$@"; do
    bash "$SCRIPTS/refine_pair.sh" "$id" || failed+=("$id")
done
echo "=== batch finished: $# pairs, ${#failed[@]} failed ${failed[*]:-}"
