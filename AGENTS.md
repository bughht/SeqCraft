# Repository agent instructions

When creating or editing anything under `examples/`, follow
[`docs/writing_examples.md`](docs/writing_examples.md).

When writing or documenting a reusable component, follow
[`docs/writing_a_module.md`](docs/writing_a_module.md) — including **Public API prose**, which
covers what belongs in a public docstring and what does not.

Examples are user-facing teaching material and docstrings are API documentation; neither is a
pull-request or design record. Development history, evidence arguments, rejected alternatives and
architecture rationale belong in `CHANGELOG.md`, the pull request, or `tools/module_mining/`.

Those two documents are the single source for those rules; this file and `CLAUDE.md` only point at
them.
