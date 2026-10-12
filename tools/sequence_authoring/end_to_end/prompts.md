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

### 1. Minimal and ambiguous

```text
Write an EPI sequence.
```

### 2. A clinical goal rather than a sequence name

```text
I need a T2-weighted 2D brain protocol at 3 T.
```

### 3. A named acquisition with structure above the repetition

```text
Build an MP2RAGE acquisition and show me the two inversion images.
```

### 4. A measured quantity rather than a picture

```text
I want to measure flow velocity in the carotid.
```

### 5. Prospective

Reserved. One genuinely prospective case is chosen later and added here **after** it has been run,
so that the first agent to see it has not also seen it sitting in the repository. Replaying an
acquisition this repository already ships is useful as a regression but does not demonstrate
independent design judgement.
