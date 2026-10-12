#!/usr/bin/env bash
#
# Build one environment for a skill end-to-end test, in which a FRESH agent session is given a
# request from `prompts.md` and nothing else.
#
#   ./setup.sh source <dir>     a SeqCraft checkout; the skill is already at .claude/skills/
#   ./setup.sh project <dir>    an independent project with SeqCraft installed and no checkout
#
# This script only builds the environment. It does not run the test, and it holds no knowledge of
# what a good answer looks like -- the point of the exercise is that the agent decides.
#
# <dir> must not exist, and must be outside the planning repository that holds the review criteria.
# SeqCraft is nested inside that repository, so a checkout placed there would let the agent reach
# them by going up one directory. Note that this is a convention, not isolation: nothing here
# restricts filesystem access, and an independent acceptance test needs an environment that does.
#
# **This script never deletes anything.** Every run wants a new directory. If setup fails partway,
# whatever it created is left in place for inspection; remove it yourself.
set -euo pipefail

MODE="${1:?usage: setup.sh source|project <dir>}"
TARGET="${2:?usage: setup.sh source|project <dir>}"

case "$MODE" in source|project) ;; *) echo "unknown mode: $MODE" >&2; exit 1 ;; esac

# ---------------------------------------------------------------- validate before touching disk
# Resolve without requiring the target to exist, and with symlinked ancestors followed, so that a
# link cannot point the "outside the planning repository" check at the wrong tree.
PARENT_RAW="$(dirname "$TARGET")"
[ -d "$PARENT_RAW" ] || { echo "refusing: parent directory does not exist: $PARENT_RAW" >&2; exit 1; }
PARENT="$(cd -P "$PARENT_RAW" && pwd -P)"
ROOT="$PARENT/$(basename "$TARGET")"

# -e is false for a dangling symlink, so test the link itself as well.
if [ -e "$ROOT" ] || [ -L "$ROOT" ]; then
  echo "refusing: $ROOT already exists; each end-to-end test uses a new directory" >&2
  exit 1
fi

# Walk every ancestor looking for the planning repository's marker directory.  Checking the
# enclosing *git* repository instead would be fooled by a nested repository, which reports its own
# top level and hides the outer one; this looks for the thing we actually care about.
PROBE="$ROOT"
while [ "$PROBE" != / ]; do
  PROBE="$(dirname "$PROBE")"
  if [ -d "$PROBE/docs/plans/sequence-authoring" ]; then
    echo "refusing: $ROOT sits under $PROBE, which holds the end-to-end review criteria" >&2
    exit 1
  fi
done

HERE="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO="${REPO:-$(cd -P "$HERE/../../.." && pwd -P)}"     # this file is in tools/sequence_authoring/end_to_end/
BRANCH="${BRANCH:-$(git -C "$REPO" rev-parse --abbrev-ref HEAD)}"

# ------------------------------------------------------------------------------------- build it
mkdir "$ROOT"                                           # not -p: a race that created it should fail

case "$MODE" in
  source)
    git clone --quiet --branch "$BRANCH" "$REPO" "$ROOT/SeqCraft"
    VENV="$ROOT/SeqCraft/.venv"
    python3 -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
    "$VENV/bin/python" -m pip install --quiet -e "$ROOT/SeqCraft[viz]"
    WHERE="$ROOT/SeqCraft"
    HOW='already present at .claude/skills/sequence-authoring/ -- native discovery is under test'
    ;;
  project)
    VENV="$ROOT/.venv"
    python3 -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
    "$VENV/bin/python" -m pip install --quiet "$REPO"
    mkdir -p "$ROOT/.claude/skills/sequence-authoring"
    cp "$REPO/.claude/skills/sequence-authoring/SKILL.md" "$ROOT/.claude/skills/sequence-authoring/"
    WHERE="$ROOT"
    HOW='copied into .claude/skills/sequence-authoring/ -- there is no packaged install for it yet'
    ;;
esac

# The agent must not be the one to discover that the environment cannot import what it needs.
if ! "$VENV/bin/python" -c 'import seqcraft, pypulseq' 2>/dev/null; then
  echo "FAILED: this environment cannot import seqcraft." >&2
  echo "  $ROOT is left in place for inspection; remove it yourself." >&2
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
