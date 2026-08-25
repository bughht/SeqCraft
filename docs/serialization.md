# Serialization — the gap, and what is deferred

**Status: the one-way MATLAB compiler-input exchange is implemented; general Python tree
serialization and `.seq` provenance sidecars remain deferred.** MATLAB can write a LogicBlock,
official `mr.opts`, official `mr.make*` events, and definitions through `seqcraft.writeLBTX`; the
Python adapter reads that file only to invoke the existing compiler. See
[ADR-005](adr/005-logicblock-tree-exchange.md) and
[`matlab_interoperability.md`](matlab_interoperability.md).

---

## What was lost

`CompiledSequence.write()` used to drop a JSON sidecar beside every `.seq`:

```text
gre_2d.seq
gre_2d.seq.json      versions, git commit + dirty flag, the definitions,
                     every field of the Opts, the achieved duration, the file's sha256
```

It went with the result wrapper. Nothing replaced it, so **a `.seq` written today does not record
what produced it.** `seq.write(path)` is pypulseq's own writer and writes the `.seq` alone.

That is a real regression and it is worth naming precisely, because the sidecar was solving a real
problem. The reference implementation's archived files carried no b-value, no diffusion directions,
no moment order, no acceleration and no partial-Fourier fraction anywhere in `[DEFINITIONS]`; those
survived only as substrings of a filename built by a forty-line `seq_file_name += ...` ladder in a
notebook cell. The parameter set that produced a given `.seq` was unrecoverable, and two files
differing only in a parameter nobody thought to put in the name were indistinguishable.

## What still works

The `[DEFINITIONS]` block does, and it is set **during the compile** rather than at write time:

```python
seq = sc.compile(tree, opts, name='gre', definitions={'FOV': [0.22, 0.22, 0.005], 'TE': 8e-3})
seq.definitions          # already populated, plus Name and TotalDuration
seq.write('gre.seq')     # nothing else has to survive until here
```

So anything you put in `definitions=` is in the file and is read back by any pulseq reader. What is
*not* there is the environment: which commit, which pypulseq, which `Opts`.

**In the meantime**, if reproducibility matters for a study, record it yourself beside the file —
`vars(opts)`, `git rev-parse HEAD`, `git status --porcelain`, and the parameters your script took.
Twenty lines, and it is the same twenty lines the sidecar was.

## Why it was not simply kept

Two reasons, and the second is the one that decided it.

**It was shaped around a type that no longer exists.** `write_sidecar` was called from
`CompiledSequence.write`, took `self.definitions`, `self.opts`, `self.n_blocks`, `self.duration_s`
and `self.report.issues`, and existed inside a `result/` package whose whole purpose was to hold
those five things. Porting it would have meant deciding, under time pressure and in the middle of a
wide refactor, what its replacement's shape is.

**The harder question was never answered.** A sidecar records what produced a `.seq`; it does not
let you *rebuild* it. The thing actually wanted is tree save/load — serialise a `LogicBlock` and
its `Opts`, reload it, and get the same sequence — and a provenance record is a weak substitute
that looks like a strong one. Designing the weak one first is how the strong one stops being
written.

## What LBTX now provides

- A versioned JSON structure for tags, ordered nodes, relative timing, scanner options, definitions,
  and official Pulseq event fields.
- MATLAB `Module` and `LogicBlock` builders using native `mr.opts` and `mr.make*` values.
- One direction only: MATLAB writes, then the Python adapter reads and compiles.

LBTX is compiler input, not a `.seq` sidecar or a mandatory Python IR. It does not currently provide
a Python writer, a MATLAB reader, a semantic hash, provenance metadata, or general round-trip
serialization. Those are deliberately deferred instead of being inferred from the existence of the
one-way adapter.

## What a `.seq` provenance replacement still has to decide

- **Where it lives.** `sc.compile` returns a `pypulseq.Sequence`, which is pypulseq's type; a
  seqcraft sidecar cannot hang off `write()` any more. Most likely a free function:
  `sc.write_with_provenance(seq, path, opts, extra=...)`, or a `provenance(opts, **extra) -> dict`
  the caller writes itself.
- **What it records.** The old set was right as far as it went. It should also carry whatever the
  *tree* knew that the `.seq` does not.
- **How it links to LBTX.** This remains undecided; the first version has no semantic hash or
  provenance contract.
- **Determinism.** The old one was careful about this and the care is worth keeping: `sort_keys`
  everywhere, no wall-clock value in the `.seq` itself, numpy arrays summarised by shape/dtype/hash
  rather than dumped, and a dirty git tree recorded as such so a comparison can be marked
  indeterminate rather than reported as a false pass.

## What is *not* being considered

A `.seq` importer. Reading a pulseq file back into a `LogicBlock` is a different problem — the file
has already had the compiler applied to it, so the tree that comes out is not the tree that went
in, and a round trip that silently changes the structure is worse than no round trip.
