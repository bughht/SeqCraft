"""
Check a candidate record against the shape the three pilots actually produced.

This is **not** a schema in the runtime sense: no classes are constructed, nothing is coerced, and
extra keys are never an error.  It reports.  The distinction is deliberate -- freezing three
candidates' vocabulary into a Pydantic model would harden the shape of exactly three examples, and
the fields that became mandatory below were found by intersecting those three rather than by
design.  A fourth candidate is expected to add keys, and the way it adds them is evidence.

What is required, and why each one survived:

    id, name, capability.{summary,family}, target.{layer,untouched}, evidence, validation,
    status, traffic_light, open_questions

    -- the thirteen field paths present in tse/, radial/ and gre3d/ candidate.yaml alike.

Four requirements are NOT from the intersection.  They are the rules the pilots learned the hard
way, and they are enforced because the pilot that needed each one did not have it:

    acceptance.does_not_establish     rule C: a claim with nothing out of scope is unscoped
    validation.layers                 a Layer-1 GREEN read as an end-to-end claim
    evidence[].independence           three references agreeing may be one witness copied twice
    reason_code                       a YELLOW nobody can triage or close
    evidence_state                    deferred work that survived only in prose, with no trigger

`evidence_state` is required only where something is actually owed -- a deferred layer, or a
non-GREEN light.  A candidate that owes nothing does not carry an empty block for symmetry.

`validation`'s interior is otherwise unconstrained on purpose: not one of its sub-keys appeared in
all three pilots, because three acceptance claims needed three different shapes of evidence.

    python tools/module_mining/schema.py path/to/candidate.yaml [more.yaml ...]
    python tools/module_mining/schema.py --inventory <records...>   # the validation-debt table

Exits non-zero if any file has an error.  Warnings never fail the run.
"""

from __future__ import annotations

import logging
import pathlib
import sys
from pathlib import Path
from typing import Any

import yaml

#: Present in tse/, radial/ and gre3d/ alike.  Dotted paths walk nested mappings.
_REQUIRED = (
    'id', 'name',
    'capability.summary', 'capability.family',
    'target.layer', 'target.untouched',
    'evidence', 'validation',
    'status', 'traffic_light', 'open_questions',
)

_FAMILIES = frozenset({'readout', 'rf', 'preparation', 'encoding', 'kernel', 'imaging'})

#: `UNDECIDED` exists because a YELLOW record still needs a status and every other value asserts
#: an action nobody took.  The historical TSE record needed it and wrote `MORE_EVIDENCE_REQUIRED`
#: by hand; a prospective run needed it again.  Pairing a provisional `NEW_LEAF` with a YELLOW
#: puts a decision in the record that was never made.
_STATUSES = frozenset({
    'NEW_LEAF', 'NEW_KERNEL', 'NEW_IMAGING', 'PROMOTE_NOTEBOOK', 'EXTEND_EXISTING',
    'NO_NEW_MODULE', 'UNDECIDED',
})

#: `APPROVED_FOR_IMPLEMENTATION` is the state between a settled design and a built module: the
#: boundary, the ownership and any mode contract have been reviewed and accepted, and nothing
#: exists yet.  GREEN cannot express it -- GREEN means the checks passed, and there is nothing to
#: check -- and YELLOW would say "stop and ask a human" when the human has already answered.
_LIGHTS = frozenset({'GREEN', 'APPROVED_FOR_IMPLEMENTATION', 'YELLOW', 'RED'})

#: reference/outcomes.md.  A YELLOW or RED outside this set cannot be triaged.
_REASON_CODES = frozenset({
    'PHYSICAL_BOUNDARY_UNCLEAR', 'MODE_CONTRACT_INCOMPLETE',
    'REFERENCES_DISAGREE', 'REFERENCE_NOT_INDEPENDENT',
    'VALIDATION_ORACLE_UNTRUSTED', 'SHARED_LEAF_IMPACT_UNRESOLVED',
    'COMPILER_CHANGE_CLAIM_NEEDS_MINIMAL_REPRODUCER',
    'WRAPPER_ONLY', 'DUPLICATES_EXISTING_MODULE', 'EXTRACTION_CHANGES_THE_PHYSICS',
})

_ROLES = frozenset({
    'primary', 'authoritative-external', 'independent-formulation', 'independent-implementation',
    'supporting', 'architecture-evidence', 'deferred',
})

_LAYER_STATES = frozenset({'GREEN', 'YELLOW', 'RED', 'DEFERRED', 'NOT_APPLICABLE'})

_DEFERRED = frozenset({'DEFERRED', 'NOT_APPLICABLE'})

#: Coarse-scan coverage, optional on a candidate record.  Separate from `status` because "what
#: already covers this" and "what we decided to do" are different questions -- "SeqCraft can build
#: it" is not "SeqCraft ships the right abstraction".  reference/outcomes.md.
_COVERAGE_CLASSES = frozenset({
    'DIRECT_SHIPPED_DUPLICATE', 'DEGENERATE_CASE_OF_EXISTING_MODULE', 'COMPOSITION_COVERED',
    'NOTEBOOK_ONLY_EXISTING', 'PRIMITIVE_COMPOSITION_COVERED',
})

#: May travel alongside a coverage class; not a verdict.
_COVERAGE_FLAGS = frozenset({'ARCHITECTURE_REVISIT_CANDIDATE'})

#: Three kinds of evidence a human supplies, deliberately not one `manual` bucket.  They differ in
#: who can do them and what they cost: the first is a workstation afternoon, the second is a
#: design decision, the third needs a magnet.
_HUMAN_STATES = frozenset({'DONE', 'NEEDED', 'NOT_NEEDED'})

#: Scanner work is claim-driven.  There is deliberately no rule that every GREEN module is
#: eventually scanner-tested -- for many leaves an independent reference comparison is more direct
#: evidence than an image, and the question is whether hardware is needed for *this* claim.
_SCANNER_STATES = frozenset({
    'NOT_REQUIRED_FOR_CLAIM', 'RECOMMENDED', 'REQUIRED', 'FUTURE', 'DONE',
})

_CONSUMER_CLASSES = frozenset({
    'NO_BEHAVIOR_CHANGE', 'NEW_EARLY_REFUSAL_FOR_PREVIOUSLY_INVALID_INPUT',
    'INTENTIONAL_BEHAVIOR_CHANGE', 'UNKNOWN / NEEDS_REVIEW',
})


def _walk(record: Any, dotted: str) -> Any:
    """Return the value at a dotted path, or ``None`` if any step is missing."""
    node = record
    for part in dotted.split('.'):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def _check_evidence(record: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    """Provenance, role vocabulary, and the independence question."""
    evidence = record.get('evidence')
    if not isinstance(evidence, list) or not evidence:
        errors.append('evidence: must be a non-empty list')
        return

    for index, item in enumerate(evidence):
        where = f'evidence[{index}]'
        if not isinstance(item, dict):
            errors.append(f'{where}: must be a mapping')
            continue
        if not item.get('repo'):
            errors.append(f'{where}: repo is required')
        if not item.get('license'):
            warnings.append(f'{where}: no license recorded')
        role = item.get('role')
        if role is None:
            errors.append(f'{where}: role is required')
        elif role not in _ROLES:
            # GRE3D wrote roles as sentences -- "primary, slab-selective -- the only witness".
            # That is a role plus a scope plus a judgement in one string; split it.
            warnings.append(
                f'{where}: role {role!r} is not in the vocabulary; if it carries scope or an '
                f'independence judgement, move those to `covers` and `independence`'
            )
        if item.get('executable') == 'executed' and not item.get('commit'):
            errors.append(f'{where}: an executed reference needs a commit')

    if len(evidence) > 1:
        # All three pilots made the independence judgement and none of them had a field for it:
        # TSE put "byte-identical modulo the import name" in a note, Radial put "ported -- not
        # independent evidence" in a note, GRE3D put it in a top-level `independence` block.  So
        # distinguish "not recorded" from "recorded somewhere unstructured" -- the second is a
        # migration, the first is missing work.
        elsewhere = isinstance(record.get('independence'), (str, dict))
        silent = [i for i, e in enumerate(evidence)
                  if isinstance(e, dict) and not e.get('independence') and not e.get('note')]
        buried = [i for i, e in enumerate(evidence)
                  if isinstance(e, dict) and not e.get('independence') and e.get('note')]
        if silent and not elsewhere:
            errors.append(
                f'evidence: {len(evidence)} references cited and '
                f'{[f"evidence[{i}]" for i in silent]} say nothing about independence -- '
                'agreement among ports of one design is one witness, not several'
            )
        elif silent:
            warnings.append(
                f'evidence: {[f"evidence[{i}]" for i in silent]} carry no per-item '
                '`independence`; a record-level `independence` block covers it less precisely'
            )
        if buried:
            warnings.append(
                f'evidence: {[f"evidence[{i}]" for i in buried]} may state independence inside '
                '`note` -- move the judgement into `independence` so it can be checked'
            )


def _check_acceptance(record: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    """Rule C, and the ordering that keeps a claim from describing the code."""
    acceptance = record.get('acceptance')
    if acceptance is None:
        warnings.append('acceptance: absent -- the three pilots recorded this in prose reports, '
                        'but a claim that is not in the record cannot be checked against')
        return
    if not acceptance.get('does_not_establish'):
        errors.append(
            'acceptance.does_not_establish: required and must be non-empty -- a comparison with '
            'nothing out of scope has not been scoped (rule C)'
        )
    if not acceptance.get('establishes'):
        errors.append('acceptance.establishes: required and must be non-empty (rule C)')
    if acceptance.get('defined_before_extraction') is not True:
        warnings.append(
            'acceptance.defined_before_extraction is not true -- an acceptance criterion written '
            'after the code exists is a description of the code'
        )


def _check_validation(record: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    """Only `layers` is constrained; the rest is family-specific by design."""
    validation = record.get('validation')
    if not isinstance(validation, dict):
        errors.append('validation: must be a mapping')
        return

    layers = validation.get('layers')
    if not isinstance(layers, dict):
        errors.append(
            'validation.layers: required -- per-layer status is what stops a Layer-1 GREEN being '
            'summarised as an end-to-end claim'
        )
        return

    for layer in ('layer_1', 'layer_2', 'layer_3'):
        state = layers.get(layer)
        if state is None:
            errors.append(f'validation.layers.{layer}: required')
        elif state not in _LAYER_STATES:
            errors.append(f'validation.layers.{layer}: {state!r} not in {sorted(_LAYER_STATES)}')

    if layers.get('layer_3') in _DEFERRED and not layers.get('layer_3_reason'):
        errors.append(
            'validation.layers.layer_3_reason: required when layer 3 is deferred -- '
            'record why the oracle is not yet trusted, not merely that it is missing'
        )

    # Rule B binds a *promoted* module.  A design that has been approved but not written has
    # nothing to emit, and demanding an inspection of it would be asking for evidence that cannot
    # exist yet.
    if record.get('traffic_light') == 'GREEN' and not validation.get('emitted_inspection'):
        modes = _walk(record, 'contract.modes')
        if modes:
            errors.append(
                f'validation.emitted_inspection: {len(modes)} physical modes declared but no '
                'emitted-sequence inspection recorded (rule B)'
            )


def _check_evidence_state(record: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    """
    The reviewer-facing block: what is owed, and what would make it actionable.

    `established` and `not_established` are NOT duplicated here -- they already live in
    `acceptance.establishes` and `acceptance.does_not_establish`, and a second copy would drift.
    What this adds is the part no field carried: each deferral's trigger, and the three kinds of
    human evidence kept apart from each other.
    """
    state = record.get('evidence_state')
    if state is None:
        light = record.get('traffic_light')
        layers = _walk(record, 'validation.layers') or {}
        owes = any(v in _DEFERRED for k, v in layers.items() if k.startswith('layer_'))
        if owes or light in {'YELLOW', 'RED'}:
            errors.append(
                'evidence_state: required -- this record defers a layer or is not GREEN, so a '
                'reviewer needs to know what evidence is owed and what would make it actionable'
            )
        return
    if not isinstance(state, dict):
        errors.append('evidence_state: must be a mapping')
        return

    for field, allowed in (('human_review', _HUMAN_STATES),
                           ('experiment_design_review', _HUMAN_STATES),
                           ('scanner_validation', _SCANNER_STATES)):
        value = state.get(field)
        if value is None:
            errors.append(f'evidence_state.{field}: required, one of {sorted(allowed)}')
        elif value not in allowed:
            errors.append(f'evidence_state.{field}: {value!r} not in {sorted(allowed)}')

    deferred = state.get('deferred')
    if deferred is None:
        errors.append('evidence_state.deferred: required; use [] when nothing is owed')
        return
    for index, item in enumerate(deferred if isinstance(deferred, list) else []):
        where = f'evidence_state.deferred[{index}]'
        if not isinstance(item, dict):
            errors.append(f'{where}: must be a mapping')
            continue
        for field in ('what', 'why', 'revisit_trigger'):
            if not item.get(field):
                errors.append(
                    f'{where}.{field}: required -- a deferral without a trigger is an omission '
                    'with better prose'
                )


def _check_dependency_impact(record: dict[str, Any], errors: list[str], warnings: list[str]) -> None:
    """Rule D fires only when the candidate changes an existing leaf."""
    extended = _walk(record, 'target.extended')
    if not extended:
        return
    impact = record.get('dependency_impact')
    if not impact:
        errors.append(
            f'dependency_impact: required -- target.extended changes {len(extended)} existing '
            'file(s), so every consumer needs classifying (rule D)'
        )
        return
    for index, entry in enumerate(impact if isinstance(impact, list) else []):
        label = (entry or {}).get('classification') if isinstance(entry, dict) else None
        if label not in _CONSUMER_CLASSES:
            errors.append(
                f'dependency_impact[{index}]: classification {label!r} not in '
                f'{sorted(_CONSUMER_CLASSES)}'
            )
        elif label == 'UNKNOWN / NEEDS_REVIEW' and record.get('traffic_light') == 'GREEN':
            errors.append(
                f'dependency_impact[{index}]: GREEN with an unresolved consumer -- '
                'this is YELLOW - SHARED_LEAF_IMPACT_UNRESOLVED'
            )


def check(record: dict[str, Any]) -> tuple[list[str], list[str]]:
    """
    Return ``(errors, warnings)`` for one parsed candidate record.

    A file carrying `evidence_state` and no `capability` is an **addendum**: the evidence state
    for a candidate whose own record is closed.  Only the addendum's block is checked, because
    the rest of the record lives -- unchanged, on purpose -- in the file it points at.
    """
    errors: list[str] = []
    warnings: list[str] = []

    if 'evidence_state' in record and 'capability' not in record:
        if not record.get('applies_to'):
            errors.append('applies_to: required on an addendum -- name the record it belongs to')
        _check_evidence_state(record, errors, warnings)
        return errors, warnings

    for dotted in _REQUIRED:
        if _walk(record, dotted) is None:
            errors.append(f'{dotted}: required (present in all three pilots)')

    family = _walk(record, 'capability.family')
    if family is not None and family not in _FAMILIES:
        errors.append(f'capability.family: {family!r} not in {sorted(_FAMILIES)}')

    status = record.get('status')
    if status is not None:
        # TSE recorded 'NEW_KERNEL + PROMOTE_NOTEBOOK'; a compound is legitimate.
        parts = {p.strip() for p in str(status).split('+')}
        unknown = parts - _STATUSES
        if unknown:
            errors.append(f'status: {sorted(unknown)} not in {sorted(_STATUSES)}')

    light = record.get('traffic_light')
    if light is not None and light not in _LIGHTS:
        errors.append(f'traffic_light: {light!r} not in {sorted(_LIGHTS)}')

    reason = record.get('reason_code')
    if light in {'YELLOW', 'RED'}:
        if not reason:
            errors.append(
                f'reason_code: required for {light} -- a traffic light without one cannot be '
                'triaged or closed (reference/outcomes.md)'
            )
        elif reason not in _REASON_CODES:
            errors.append(f'reason_code: {reason!r} not in the catalogue')
    elif light == 'GREEN' and reason:
        warnings.append(f'reason_code: {reason!r} set on a GREEN record')

    if status == 'UNDECIDED' and light in {'GREEN', 'APPROVED_FOR_IMPLEMENTATION'}:
        errors.append(f'status: UNDECIDED cannot be {light} -- both are decisions')
    if light == 'APPROVED_FOR_IMPLEMENTATION':
        layers = _walk(record, 'validation.layers') or {}
        built = [k for k, v in layers.items() if k.startswith('layer_') and v == 'GREEN']
        if built:
            warnings.append(
                f'traffic_light: APPROVED_FOR_IMPLEMENTATION with {sorted(built)} already GREEN '
                '-- if the module exists and its layers pass, this is GREEN'
            )

    coverage = record.get('coverage')
    if coverage is not None:
        parts = {c.strip() for c in str(coverage).split('+')}
        unknown = parts - _COVERAGE_CLASSES - _COVERAGE_FLAGS
        if unknown:
            errors.append(f'coverage: {sorted(unknown)} not in {sorted(_COVERAGE_CLASSES)}')
        elif not parts & _COVERAGE_CLASSES:
            errors.append('coverage: a flag alone is not a classification')

    if record.get('open_questions') == []:
        warnings.append('open_questions: empty -- all three pilots kept some, even at GREEN')

    if _walk(record, 'contract') is None and _walk(record, 'physics') is None:
        warnings.append(
            'contract: absent -- TSE and Radial called this `physics` and GRE3D called it '
            '`working_hypothesis`; the skill names it `contract`'
        )

    _check_evidence(record, errors, warnings)
    _check_acceptance(record, errors, warnings)
    _check_validation(record, errors, warnings)
    _check_evidence_state(record, errors, warnings)
    _check_dependency_impact(record, errors, warnings)
    return errors, warnings


def _inventory(paths: list[Path]) -> None:
    """
    Print the cross-candidate validation-debt table, read from the records.

    Generated rather than maintained: a hand-written inventory beside the records is a second
    copy of the same facts, and the copy is the one that goes stale.
    """
    rows: list[tuple[str, str, str, str, int]] = []
    addenda: dict[str, dict[str, Any]] = {}
    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        try:
            record = yaml.safe_load(path.read_text())
        except (OSError, yaml.YAMLError):
            continue
        if not isinstance(record, dict):
            continue
        if 'evidence_state' in record and 'capability' not in record:
            addenda[str(record.get('applies_to', ''))] = record
        else:
            records[str(path)] = record

    for path, record in sorted(records.items()):
        state = record.get('evidence_state')
        for target, addendum in addenda.items():
            if target and pathlib.Path(target).name == pathlib.Path(path).name \
                    and pathlib.Path(target).parent.name == pathlib.Path(path).parent.name:
                state = addendum['evidence_state']
        state = state or {}
        deferred = state.get('deferred') or []
        name = str(record.get('name') or '')
        # A candidate whose working name is still open says so in `name`; the id is the label.
        if name in {'', 'UNDECIDED'}:
            name = str(record.get('id', '?'))
        rows.append((
            name,
            f"{record.get('status', '?')} / {record.get('traffic_light', '?')}",
            str(state.get('scanner_validation', '--')),
            str(state.get('human_review', '--')),
            len(deferred),
        ))

    width = max((len(r[0]) for r in rows), default=4)
    logging.info('')
    logging.info('%-*s  %-28s  %-24s  %-10s  %s', width, 'candidate', 'status / light',
                 'scanner', 'human', 'owed')
    for name, status, scanner, human, owed in rows:
        logging.info('%-*s  %-28s  %-24s  %-10s  %d', width, name, status, scanner, human, owed)
    logging.info('')
    for path, record in sorted(records.items()):
        state = record.get('evidence_state') or {}
        for target, addendum in addenda.items():
            if target and pathlib.Path(target).name == pathlib.Path(path).name \
                    and pathlib.Path(target).parent.name == pathlib.Path(path).parent.name:
                state = addendum['evidence_state']
        label = str(record.get('name') or '')
        if label in {'', 'UNDECIDED'}:
            label = str(record.get('id', '?'))
        for item in state.get('deferred') or []:
            logging.info('%s owes: %s', label,
                         ' '.join(str(item.get('what', '')).split())[:96])
            logging.info('   trigger: %s',
                         ' '.join(str(item.get('revisit_trigger', '')).split())[:96])


def main(argv: list[str]) -> int:
    """Validate each file named on the command line; return 1 if any had an error."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    if not argv:
        logging.info('%s', __doc__)
        return 2

    if argv and argv[0] == '--inventory':
        _inventory([Path(name) for name in argv[1:]])
        return 0

    worst = 0
    for name in argv:
        path = Path(name)
        try:
            record = yaml.safe_load(path.read_text())
        except (OSError, yaml.YAMLError) as error:
            logging.info('%s: cannot read -- %s', path, error)
            worst = 1
            continue
        if not isinstance(record, dict):
            logging.info('%s: top level is not a mapping', path)
            worst = 1
            continue

        errors, warnings = check(record)
        logging.info('')
        logging.info('%s', path)
        for message in errors:
            logging.info('  ERROR    %s', message)
        for message in warnings:
            logging.info('  warning  %s', message)
        if not errors and not warnings:
            logging.info('  ok')
        elif not errors:
            logging.info('  ok, %d warning(s)', len(warnings))
        worst = max(worst, 1 if errors else 0)
    return worst


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
