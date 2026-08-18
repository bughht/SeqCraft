"""Magnetisation preparation: the folder rule is ``rf.use`` in {inversion, saturation, preparation}.

``rf/`` keeps ``excitation`` and ``refocusing``, which are the two uses that belong to the imaging
train itself; everything played *before* it is here.  Classes are named ``<physics>Prep``, because
several distinct physics share one ``use`` and the role alone would not tell them apart -- see
``ir_prep.py``.

The folder is taxonomy only: ``seqcraft.modules``
imports flat, so nothing outside the package writes this path.
"""
