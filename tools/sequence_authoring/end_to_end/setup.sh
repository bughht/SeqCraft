#!/usr/bin/env bash
#
# Build one environment for a skill end-to-end test, in which a FRESH agent session is given a
# request from `prompts.md` and nothing else.
#
#   ./setup.sh source <dir>     a SeqCraft checkout; the skill is already at .claude/skills/
#   ./setup.sh project <dir>    an independent project with SeqCraft installed and no checkout
#
# This script only builds the environment. It does not run the test, and it must not acquire any
# knowledge of what a good answer looks like -- the point of the exercise is that the agent decides.
#
# <dir> must be outside the planning repository that holds the review criteria. SeqCraft is nested
# inside it, so a checkout placed there would let the agent reach them by going up one directory.
# Note that this is a convention, not isolation: nothing here restricts filesystem access, and an
# independent acceptance test needs an execution environment that does.
set -euo pipefail

MODE="${1:?usage: setup.sh source|project <dir>}"
TARGET="${2:?usage: setup.sh source|project <dir>}"
ROOT="$(cd "$(dirname "$TARGET")" && pwd)/$(basename "$TARGET")"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${REPO:-$(cd "$HERE/../../.." && pwd)}"          # this file lives in tools/sequence_authoring/end_to_end/
BRANCH="${BRANCH:-$(git -C "$REPO" rev-parse --abbrev-ref HEAD)}"
EXTRAS="${EXTRAS:-}"                                   # e.g. EXTRAS='[rf]' before PR #54 merges

# Refuse a target inside any repository that keeps the review criteria beside it.  The check walks
# up to the nearest ancestor that exists: $ROOT itself usually does not yet, and an earlier version
# asked git about a path containing a missing component, which failed silently and let the guard
# through -- it then installed 224 MB into the repository it was meant to protect.
ANCESTOR="$ROOT"
while [ ! -d "$ANCESTOR" ] && [ "$ANCESTOR" != / ]; do ANCESTOR="$(dirname "$ANCESTOR")"; done
if OUTER="$(git -C "$ANCESTOR" rev-parse --show-toplevel 2>/dev/null)"; then
  if [ -d "$OUTER/docs/plans/sequence-authoring" ]; then
    echo "refusing: $ROOT sits inside $OUTER, which holds the end-to-end review criteria" >&2
    exit 1
  fi
fi

rm -rf "$ROOT"
mkdir -p "$ROOT"

case "$MODE" in
  source)
    git clone --quiet --branch "$BRANCH" "$REPO" "$ROOT/SeqCraft"
    VENV="$ROOT/SeqCraft/.venv"
    python3 -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
    "$VENV/bin/python" -m pip install --quiet -e "$ROOT/SeqCraft[rf,viz]"
    WHERE="$ROOT/SeqCraft"
    HOW='already present at .claude/skills/sequence-authoring/ -- native discovery is under test'
    ;;
  project)
    VENV="$ROOT/.venv"
    python3 -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
    "$VENV/bin/python" -m pip install --quiet "${REPO}${EXTRAS}"
    mkdir -p "$ROOT/.claude/skills/sequence-authoring"
    cp "$REPO/.claude/skills/sequence-authoring/SKILL.md" "$ROOT/.claude/skills/sequence-authoring/"
    WHERE="$ROOT"
    HOW='copied into .claude/skills/sequence-authoring/ -- there is no packaged install for it yet'
    ;;
  *) echo "unknown mode: $MODE" >&2; exit 1 ;;
esac

# The agent must not be the one to discover that the environment cannot import what it needs.
if ! "$VENV/bin/python" -c 'import seqcraft, pypulseq' 2>/dev/null; then
  echo "FAILED: this environment cannot import seqcraft." >&2
  echo "  A base install could not, until SigPy became a base dependency; set EXTRAS='[rf]'." >&2
  exit 1
fi
"$VENV/bin/python" - <<'CHECK'
import importlib.metadata as md
import seqcraft as sc
mods = [n for n in dir(sc.modules) if not n.startswith('_')]
print(f"  seqcraft {sc.__version__} | pypulseq {md.version('pypulseq')} | {len(mods)} modules")
CHECK

echo "  work in:  $WHERE"
echo "  skill:    $HOW"
echo "  prompt:   give the agent one line from prompts.md, verbatim, and nothing else"
