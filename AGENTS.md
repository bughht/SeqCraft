# Repository agent instructions

When creating or editing anything under `examples/`, follow
[`docs/writing_examples.md`](docs/writing_examples.md).

When writing or documenting a reusable component, follow
[`docs/writing_a_module.md`](docs/writing_a_module.md) — including **Public API prose**, which
covers what belongs in a public docstring and what does not.

When writing a complete MRI acquisition with SeqCraft — choosing a protocol, building it,
validating the emitted `.seq` and reporting its limits — follow
[`.claude/skills/sequence-authoring/SKILL.md`](.claude/skills/sequence-authoring/SKILL.md). That
path is Claude Code's skill-discovery location; this line is a reference to the same single copy,
not a second discovery mechanism, and the cross-agent arrangement is still to be settled.

Examples are user-facing teaching material and docstrings are API documentation; neither is a
pull-request or design record. Development history, evidence arguments, rejected alternatives and
architecture rationale belong in `CHANGELOG.md`, the pull request, or `tools/module_mining/`.

Each of those documents is the single source for its own rules; this file and `CLAUDE.md` only
point at them.
