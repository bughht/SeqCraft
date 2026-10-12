#!/usr/bin/env bash
#
# Build one environment for a skill end-to-end test, in which a FRESH agent session is given a
# request from `prompts.md` and nothing else.
#
#   ./setup.sh source <dir>     a SeqCraft checkout; the skill is already at .claude/skills/
#   ./setup.sh project <dir>    an independent project with SeqCraft installed and no checkout
#
# By default each run builds its own virtual environment, which is the clean-install option and the
# one to use when dependency isolation is itself what you are checking. Set CONDA_ENV to reuse a
# dedicated conda environment instead, so four test cases do not rebuild four copies of torch:
#
#   CONDA_ENV=seqcraft-e2e ./setup.sh project ~/skill-tests/epi-01
#
# The working directory is always new either way -- only the interpreter is shared.
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
CONDA_ENV="${CONDA_ENV:-}"
REINSTALL="${REINSTALL:-}"                              # set to force a reinstall into a conda env

# Resolve the interpreter before anything is built, so a missing environment fails here.
if [ -n "$CONDA_ENV" ]; then
  command -v conda >/dev/null || { echo "refusing: CONDA_ENV set but conda is not on PATH" >&2; exit 1; }
  CONDA_PREFIX_DIR="$(conda env list | awk -v n="$CONDA_ENV" '$1==n {print $NF}')"
  [ -n "$CONDA_PREFIX_DIR" ] || {
    echo "refusing: no conda environment named $CONDA_ENV" >&2
    echo "  create it first, e.g.: conda create -n $CONDA_ENV python=3.12" >&2
    exit 1; }
  PY="$CONDA_PREFIX_DIR/bin/python"
  [ -x "$PY" ] || { echo "refusing: no interpreter at $PY" >&2; exit 1; }
fi

# ------------------------------------------------------------------------------------- build it
mkdir "$ROOT"                                           # not -p: a race that created it should fail

case "$MODE" in
  source)
    git clone --quiet --branch "$BRANCH" "$REPO" "$ROOT/SeqCraft"
    if [ -z "$CONDA_ENV" ]; then
      PY="$ROOT/SeqCraft/.venv/bin/python"
      python3 -m venv "$ROOT/SeqCraft/.venv"
      "$PY" -m pip install --quiet --upgrade pip
      "$PY" -m pip install --quiet -e "$ROOT/SeqCraft[viz]"
    fi
    WHERE="$ROOT/SeqCraft"
    HOW='already present at .claude/skills/sequence-authoring/ -- native discovery is under test'
    ;;
  project)
    # An installed distribution, never an editable install pointing at the development checkout --
    # the point of this mode is that no source tree is reachable from the project.
    if [ -z "$CONDA_ENV" ]; then
      PY="$ROOT/.venv/bin/python"
      python3 -m venv "$ROOT/.venv"
      "$PY" -m pip install --quiet --upgrade pip
      "$PY" -m pip install --quiet "$REPO"
    elif [ -n "$REINSTALL" ] || ! "$PY" -c 'import seqcraft' 2>/dev/null; then
      "$PY" -m pip install --quiet --upgrade "$REPO"
    fi
    mkdir -p "$ROOT/.claude/skills/sequence-authoring"
    cp "$REPO/.claude/skills/sequence-authoring/SKILL.md" "$ROOT/.claude/skills/sequence-authoring/"
    WHERE="$ROOT"
    HOW='copied into .claude/skills/sequence-authoring/ -- there is no packaged install for it yet'
    ;;
esac

# The agent must not be the one to discover that the environment cannot import what it needs.
if ! "$PY" -c 'import seqcraft, pypulseq' 2>/dev/null; then
  echo "FAILED: this environment cannot import seqcraft." >&2
  echo "  $ROOT is left in place for inspection; remove it yourself." >&2
  exit 1
fi

# Record what was actually resolved, not what was intended.  An editable install pointing back at
# the development checkout would defeat the project mode, and the only way to know is to look.
"$PY" - "$MODE" "$REPO" > "$ROOT/environment.txt" <<'CHECK'
import importlib.metadata as md
import pathlib
import sys

import pypulseq
import seqcraft as sc

mode, repo = sys.argv[1], pathlib.Path(sys.argv[2]).resolve()
where = pathlib.Path(sc.__file__).resolve().parent
editable = repo in where.parents or where == repo / 'src' / 'seqcraft'

print(f'python            {sys.executable}')
print(f'version           {sys.version.split()[0]}')
print(f'seqcraft          {sc.__version__}')
print(f'seqcraft at       {where}')
print(f'editable install  {"yes" if editable else "no"}')
print(f'pypulseq          {md.version("pypulseq")}')
print(f'pypulseq at       {pathlib.Path(pypulseq.__file__).resolve().parent}')
print(f'numpy             {md.version("numpy")}')
print(f'sigpy             {md.version("sigpy")}')
print(f'modules           {len([n for n in dir(sc.modules) if not n.startswith("_")])}')
for extra in ('MRzeroCore', 'torch', 'matplotlib'):
    try:
        print(f'{extra:<18}{md.version(extra)}')
    except md.PackageNotFoundError:
        print(f'{extra:<18}not installed')
if mode == 'project' and editable:
    raise SystemExit('FAILED: project mode resolved to an editable install of the checkout')
CHECK
sed 's/^/  /' "$ROOT/environment.txt"

echo "  ---"
echo "  work in:  $WHERE"
echo "  skill:    $HOW"
echo "  recorded: $ROOT/environment.txt"
echo "  prompt:   give the agent one line from prompts.md, verbatim, and nothing else"
