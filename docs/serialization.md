# Serialization — the gap, and what is deferred

**Status: tree exchange is implemented; `.seq` provenance sidecars remain deferred.** LBTX 0.1 can
save and load a `LogicBlock`, its scanner options, and definitions through `sc.write_logicblock`
and `sc.read_logicblock`. See [ADR-005](adr/005-logicblock-tree-exchange.md) and
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

- Versioned JSON Schema and explicit SI event/scanner fields.
- Python semantic round-trip with insertion order, nested tags, relative timing, and overlap intact.
- Separate provenance metadata and a semantic hash that excludes provenance and extensions.
- MATLAB read/write/build support and a structured Python compilation CLI.

LBTX is compiler input, not a `.seq` sidecar. It lets a caller rebuild the sequence, but it does not
automatically record the environment whenever `pypulseq.Sequence.write()` is called.

## What a `.seq` provenance replacement still has to decide

- **Where it lives.** `sc.compile` returns a `pypulseq.Sequence`, which is pypulseq's type; a
  seqcraft sidecar cannot hang off `write()` any more. Most likely a free function:
  `sc.write_with_provenance(seq, path, opts, extra=...)`, or a `provenance(opts, **extra) -> dict`
  the caller writes itself.
- **What it records.** The old set was right as far as it went. It should also carry whatever the
  *tree* knew that the `.seq` does not.
- **How it links to LBTX.** A sidecar can carry the LBTX semantic hash and producer environment
  without embedding a second copy of the tree.
- **Determinism.** The old one was careful about this and the care is worth keeping: `sort_keys`
  everywhere, no wall-clock value in the `.seq` itself, numpy arrays summarised by shape/dtype/hash
  rather than dumped, and a dirty git tree recorded as such so a comparison can be marked
  indeterminate rather than reported as a false pass.

## What is *not* being considered

A `.seq` importer. Reading a pulseq file back into a `LogicBlock` is a different problem — the file
has already had the compiler applied to it, so the tree that comes out is not the tree that went
in, and a round trip that silently changes the structure is worse than no round trip.
