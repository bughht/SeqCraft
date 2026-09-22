# Domain evidence cards

Curated notes from **authoritative domain literature** — handbooks, major reviews, classic papers
— for candidates that reach a fine scan.

They exist because the scan corpus in `../sources.yaml` is code, and code is the wrong evidence
class for the question *"is this a mature, well-established physical abstraction?"*. A family can
be textbook physics and still appear in only one of nine registered repositories.

```text
domain reference   handbook / authoritative review / classic paper
                   -> canonical concepts, terminology, semantic physical quantities,
                      and mature family boundaries

design witness     an external implementation  (../sources.yaml)
                   -> corroborates a particular realisation, timing convention or mode

oracle             analytic calculation / Bloch simulation / independent measurement
                   -> tests whether OUR realisation produces the claimed physics
```

**The number of executable repositories is not a vote on whether a physical abstraction exists.**
Code witnesses constrain implementation-specific claims; they do not decide the family.

## Rules

- One card per candidate, **one to two pages**. Not a literature review, not a database.
- Every substantive statement carries provenance: source, chapter/section, page.
- Every statement says **what it establishes** and **what it does not**.
- Sources are **not** committed. Cards cite; they do not reproduce.
- Escalate progressively: handbook or major review first; a classic original paper only if that
  leaves a question open; application-specific papers only for a concrete unresolved question.
- If a paper is relevant but the full text is unavailable, say so — what was available, what can
  safely be concluded, and what still needs the full text. **Do not reconstruct a source from
  general knowledge.**

## Cards

| card | candidate | primary source |
|---|---|---|
| [`t2prep.md`](t2prep.md) | T2 preparation (Phase C, first fine scan) | Bernstein, King & Zhou, *Handbook of MRI Pulse Sequences* (2004), §17.4 |
| [`bssfp.md`](bssfp.md) | balanced SSFP (Phase C, second fine scan) | the same, §14.1 and Table 14.2 |
