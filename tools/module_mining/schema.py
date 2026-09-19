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

`validation`'s interior is otherwise unconstrained on purpose: not one of its sub-keys appeared in
all three pilots, because three acceptance claims needed three different shapes of evidence.

    python tools/module_mining/schema.py path/to/candidate.yaml [more.yaml ...]

Exits non-zero if any file has an error.  Warnings never fail the run.
"""

from __future__ import annotations

import logging
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

_STATUSES = frozenset({
    'NEW_LEAF', 'NEW_KERNEL', 'NEW_IMAGING', 'PROMOTE_NOTEBOOK', 'EXTEND_EXISTING',
    'NO_NEW_MODULE',
})

_LIGHTS = frozenset({'GREEN', 'YELLOW', 'RED'})

#: reference/outcomes.md.  A YELLOW or RED outside this set cannot be triaged.
_REASON_CODES = frozenset({
    'PHYSICAL_BOUNDARY_UNCLEAR', 'REFERENCES_DISAGREE', 'REFERENCE_NOT_INDEPENDENT',
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

    if record.get('traffic_light') == 'GREEN' and not validation.get('emitted_inspection'):
        modes = _walk(record, 'contract.modes')
        if modes:
            errors.append(
                f'validation.emitted_inspection: {len(modes)} physical modes declared but no '
                'emitted-sequence inspection recorded (rule B)'
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
    """Return ``(errors, warnings)`` for one parsed candidate record."""
    errors: list[str] = []
    warnings: list[str] = []

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
    _check_dependency_impact(record, errors, warnings)
    return errors, warnings


def main(argv: list[str]) -> int:
    """Validate each file named on the command line; return 1 if any had an error."""
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    if not argv:
        logging.info('%s', __doc__)
        return 2

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
