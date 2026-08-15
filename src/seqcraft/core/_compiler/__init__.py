"""Private contracts and stages used by :mod:`seqcraft.core.compiler`.

This package is an implementation boundary, not a public API.  The compatible compiler entry
point remains :func:`seqcraft.core.compiler.compile_sequence` while the refactor moves one stage
at a time behind that facade.
"""

from __future__ import annotations
