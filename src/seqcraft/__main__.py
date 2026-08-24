"""Allow ``python -m seqcraft`` to use the installed command-line interface."""

from .cli import main

raise SystemExit(main())
