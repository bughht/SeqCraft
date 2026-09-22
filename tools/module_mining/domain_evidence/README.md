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

## What a card contains, and what it must not

A card is **evidence**, not environment state. Two developers with the same book on different
machines must be able to use the same card unchanged.

```text
keep                                    do not include
----                                    --------------
the bibliographic source                any absolute or local filesystem path
chapter / section                       where a PDF happens to be stored
book and page numbers                   whether a local copy is committed
the exact scope actually read           machine-specific access instructions
claims established / not established
source-native terminology
a stable source ID
```

The "only these pages were read" notes **stay** — those describe evidence scope, not file location.

## Getting at a source on your machine

Cards name a source by a stable ID. Resolution is machine-local and lives outside them:

```text
tools/module_mining/local_sources.example.json   tracked -- the template
tools/module_mining/local_sources.json           gitignored -- your paths
```

The workflow and the filename are shared; the resolved path is not, which is why exactly one
`.gitignore` entry exists for it. Resolution order:

```text
1. the source's env_override, e.g. SEQCRAFT_HANDBOOK_PDF   -- highest priority
2. a valid path in local_sources.json
3. a path the user has just supplied -> write it into local_sources.json from the template
4. otherwise: report that the source is unavailable and say which one is needed
```

If a configured path no longer exists, **report the stale mapping.** Do not search the filesystem
for a replacement.

## Rules

- One card per candidate, **one to two pages**. Not a literature review, not a database.
- Every substantive statement carries provenance: source, chapter/section, page.
- Every statement says **what it establishes** and **what it does not**.
- Cards cite; they do not reproduce. Source files live wherever the reader keeps them.
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
| [`gradient_moments.md`](gradient_moments.md) | flow encoding + moment nulling (Phase C, third fine scan) | the same, §9.2 and §10.4 |
