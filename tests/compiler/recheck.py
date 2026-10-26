"""
Re-run a real compile's against-tree check with one input changed.

A tampered expectation stands in for a compiler that made the opposite error: checking a clean
sequence against a tree that also asks for a sliver is what a compile that dropped the sliver
would look like.  To mean that, the check has to be the one ``sc.compile`` runs, on the inputs it
ran with -- the ready blocks measured before emission, the built sequence, the reported resamples
-- so these helpers capture them from a real compile instead of assembling their own.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import seqcraft as sc
import seqcraft.compiler as compiler_package
from seqcraft.compiler.verification import verify_against_tree

if TYPE_CHECKING:
    import pypulseq as pp

__all__ = ['compile_capturing', 'recheck']


def compile_capturing(
    tree: sc.LogicBlock, opts: pp.Opts, *, check: bool = True, **kwargs: Any,
) -> tuple[pp.Sequence, dict[str, Any]]:
    """
    Compile `tree`, keeping everything ``sc.compile`` handed ``verify_against_tree``.

    With ``check=False`` the check is captured but not run, to measure what it measures.
    """
    seen: dict[str, Any] = {}
    original = compiler_package.verify_against_tree

    def capture(placed: Any, targets: Any, **inputs: Any) -> None:
        seen.update(placed=placed, targets=targets, **inputs)
        if check:
            original(placed, targets, **inputs)

    compiler_package.verify_against_tree = capture
    try:
        seq = sc.compile(tree, opts, **kwargs)
    finally:
        compiler_package.verify_against_tree = original
    return seq, seen


def recheck(seen: dict[str, Any], **changes: Any) -> None:
    """Run the check again on the captured inputs, with `changes` replacing some of them."""
    inputs = {**seen, **changes}
    verify_against_tree(inputs.pop('placed'), inputs.pop('targets'), **inputs)
