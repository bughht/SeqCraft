"""
``sc.modules`` -- the concrete MR building blocks, re-exported flat.

seqcraft's core ships no recipes: the tree, the compiler and the contract are the package, and a
recipe is somebody else's sequence choices baked into library code.  This is the other half of
that argument.  A *module* is not a recipe -- it is one piece of MR physics with its arithmetic
attached, and the difference is whether the caller keeps the sequence choices.  ``GRE2D`` takes
the list of lines to acquire rather than an acceleration factor, so the sampling pattern stays
the caller's; ``CartesianLine`` computes a prephaser that cancels the ramp, which is arithmetic
nobody should have to get right twice.

Nothing here was designed in the abstract.  Every module was **extracted from working code**:
written as raw pypulseq in ``examples/gre_2d/01_build.ipynb``, simulated until the image was
right, and only then lifted out with the compiled output held fixed.  A module that cannot be
extracted without altering the sequence is not a module; one whose extraction does not shorten
the notebook is a wrapper.

The folders, and what each one means
------------------------------------
::

    modules/
      spoiler.py           spoiler()          a function, not a Module
      rf/          excitation.py       Excitation
                   refocusing.py       Refocusing
      preparation/ ir_prep.py          IRPrep
      encoding/    phase_encoding.py   PhaseEncode
      readout/     cartesian_line.py   CartesianLine
      kernel/      gre_2d_tr.py        GRE2DTR      composes leaves; one repeating unit
      imaging/     gre_2d.py           GRE2D        composes kernels; a complete scan

+-----------------+----------------------------------------------------------------+
| ``rf/``         | ``rf.use`` in {excitation, refocusing}                          |
+-----------------+----------------------------------------------------------------+
| ``preparation/``| ``rf.use`` in {inversion, saturation, preparation} -- played     |
|                 | before the imaging train                                        |
+-----------------+----------------------------------------------------------------+
| ``encoding/``   | gradients, no ADC, imposing a phase you intend to sample        |
+-----------------+----------------------------------------------------------------+
| ``readout/``    | contains an ADC                                                 |
+-----------------+----------------------------------------------------------------+
| ``kernel/``     | composes modules from more than one leaf folder -- **the        |
|                 | repeating unit**                                                |
+-----------------+----------------------------------------------------------------+
| ``imaging/``    | composes kernels -- **a complete scan**                         |
+-----------------+----------------------------------------------------------------+
| top level       | what is not an ``sc.Module`` subclass                           |
+-----------------+----------------------------------------------------------------+

The leaf rules discriminate leaf modules from one another; they cannot discriminate a leaf from a
composite, because a GRE kernel contains an ADC and so would a whole DTI scan.  The last two
rules do, and the split matches how the pieces are actually reused.

``rf/`` and ``preparation/`` split one field rather than two ideas, and the split is worth having
because the names differ in kind: ``rf.use`` **is** the role for an excitation or a refocusing, so
those classes need nothing added, while several distinct physics share ``'preparation'`` and one
shares ``'inversion'``, so those carry a ``Prep`` suffix to say which.  ``ir_prep.py`` argues it.

**Folders never appear in an import path.**  This module re-exports flat, so the taxonomy stays
cheap to revise -- which matters, because §5 of the plan that built it schedules a decision about
whether five folders still fit once there are eight or nine modules.

``sc.modules`` is not lazily resolved.  It imports pypulseq and :mod:`seqcraft.design`, both of
which ``import seqcraft`` has already paid for, and nothing else.

Examples
--------
>>> import pypulseq as pp
>>> import seqcraft as sc
>>> opts = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
...                rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
>>> gre = sc.modules.GRE2D(opts=opts, fov_mm=(250.0, 250.0), matrix=(64, 32),
...                        thickness_mm=5.0, flip_deg=15.0)
>>> seq = sc.compile(gre(lines=range(32)), opts, name='gre_2d')
>>> len(seq.block_events), round(seq.duration()[0] * 1e3, 2)
(128, 306.24)

The provenance path is the tree, and not one tag string was written:

>>> sorted({path for _, _, path in sc.flatten(gre(lines=[16]))})
[('GRE2D', 'GRE2DTR'), ('GRE2D', 'GRE2DTR', 'CartesianLine'), ('GRE2D', 'GRE2DTR', 'Excitation'), ('GRE2D', 'GRE2DTR', 'PhaseEncode'), ('GRE2D', 'GRE2DTR', 'spoiler')]
"""

from __future__ import annotations

from .encoding.phase_encoding import PhaseEncode
from .imaging.gre_2d import GRE2D
from .kernel.gre_2d_tr import GRE2DTR
from .preparation.ir_prep import IRPrep
from .readout.cartesian_line import CartesianLine
from .rf.excitation import Excitation
from .rf.refocusing import Refocusing
from .spoiler import spoiler

__all__ = [
    'CartesianLine', 'Excitation', 'GRE2D', 'GRE2DTR', 'IRPrep', 'PhaseEncode', 'Refocusing',
    'spoiler',
]
