# Sequence-authoring skill end-to-end test prompts

**Inputs only.** These are the requests used in end-to-end tests of the `sequence-authoring`
skill. The file deliberately contains no expected answers, no scoring and no verdicts: the agent
under test can read this repository, so anything written here is part of its context.

Review criteria and recorded findings live in the workspace planning repository, under
`docs/plans/sequence-authoring/`.

## What kind of test this is

```text
unit / integration tests   deterministic Python tests of SeqCraft, emitted .seq artifacts and
                           physical invariants.  tests/sequence_authoring/ holds the ones that
                           belong to this skill, and pytest runs them.

skill end-to-end tests     a FRESH agent session is given one of the requests below, uses the
                           skill, writes a complete acquisition, validates it, and reports.
                           A person reviews what comes back.  This file feeds those.

cross-agent tests          the same end-to-end test run separately under Claude Code and Codex,
                           against one copy of the skill.

benchmark                  a future controlled comparison of agent configurations.  Not part
                           of v0.
```

Nothing in `tests/` runs the prompts below, and a passing pytest run is not evidence about any of
them.

## Two things an end-to-end run measures

Keep them apart when reporting, because they fail for different reasons and have different fixes:

```text
skill discovery        does the agent find and invoke the skill from an ordinary request?
skill-guided authoring once the skill is loaded, is the acquisition better?
```

A run where the agent never loaded the skill says nothing about the second.

## How to run one

Start a **new session**. Give the agent the prompt verbatim and nothing else — do not name the
acquisition family, supply parameters it did not ask for, or point it at an existing example. The
point is to see what it does with an underspecified request, including whether and when it asks.

Let it write its artifacts into its own project directory.

## The prompts

Four cases, one per authoring capability. They are given one at a time, each to its own fresh
session, and the first is inspected before the rest run so that one environment problem does not
spoil four results.

Reusing a shipped Module correctly is a success, not a shortcut. An agent that reimplements
something the package already provides has done worse, not better, even if its version works.
Nothing below names a Module, because discovering the right abstraction is part of what is being
observed.

### A — reuse

```text
Write an EPI sequence.
```

### B — extend

```text
Write an EPI sequence that acquires the entire k-space trajectory as one continuous readout.
```

### C — compose

```text
Write a 2D T2-prepared gradient-echo acquisition.
```

### D — design

```text
Design a 2D gradient-echo acquisition that samples k-space along a rosette trajectory.
```

### Later, prospective

Reserved. One genuinely prospective case is chosen later and added **after** it has been run, so
the first agent to see it has not also seen it sitting in the repository. Replaying an acquisition
this repository already ships is useful as a regression but does not demonstrate independent design
judgement.
