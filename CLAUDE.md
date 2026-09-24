# Working in this repository

## Examples

When creating or editing anything under `examples/`, follow the notebook-writing rules in
[`docs/writing_examples.md`](docs/writing_examples.md).

## Modules and public docstrings

When writing or documenting a reusable component, follow
[`docs/writing_a_module.md`](docs/writing_a_module.md) — including **Public API prose**, which
covers what a public docstring documents and what belongs elsewhere.

## Where the rest goes

Examples are user-facing teaching material and docstrings are API documentation; neither is a
pull-request or design record. Development history, evidence arguments, rejected alternatives and
architecture rationale belong in `CHANGELOG.md`, the pull request, or
`tools/module_mining/plans/`.

Those two documents are the single source for those rules. This file, `AGENTS.md` and the
module-mining skill only point at them.
