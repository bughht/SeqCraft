# Sequence-authoring calibration prompts

**Inputs only.** This file is the set of requests used to calibrate the `sequence-authoring` skill.
It deliberately contains no expected answers, no scoring, and no verdicts: the agent being
evaluated can read this repository, so anything written here is part of its context.

The review criteria and the recorded outcomes live in the workspace planning repository, under
`docs/plans/sequence-authoring/`, and are not mirrored here.

**These are not tests.** Nothing in `tests/` runs them. They are run by hand against an agent, and
what comes back is judged by a person. The deterministic checks on an emitted acquisition live in
[`tests/sequence_authoring/`](../../../tests/sequence_authoring/) and answer a different question —
whether a file is the acquisition it claims to be, not whether the acquisition was well chosen.

## How to run one

Give the agent the prompt verbatim and nothing else. In particular do not name the acquisition
family, supply parameters it did not ask for, or point it at an existing example — the point is to
see what it does with an underspecified request, including whether and when it asks.

Let it write its artifacts under `scratch/sequence_authoring/<case>/`, which is git-ignored.

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

Reserved. One genuinely prospective case is chosen at calibration time and added here **after** it
has been run, so that the first agent to see it has not also seen it sitting in the repository. A
replay of an acquisition this repository already ships is useful as a regression but does not
demonstrate independent design judgement.
