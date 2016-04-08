DYNAMO - Quantum Dynamic Optimization Package

Version 1.4.0
Released 2015-12-??


Introduction
============

DYNAMO is a flexible framework for solving numerical optimal control
problems for quantum systems. It requires MATLAB R2011a or newer (matlab.mixin.Copyable).

The latest development version can be downloaded from our Git
repository at https://github.com/smite/Dynamo
The user manual can be found in docs/dynamo_manual.tex, together with
the bibliography.

If you use DYNAMO in your research, please add an attribution in the form of the following reference:
[TODO: DYNAMO 2 paper]


Getting started
===============

Initialize the package by running init.m

The best way to learn how to use DYNAMO is to review the demos in the
examples/ directory.


Design
======

The design is modularized and easily extendable.

DYNAMO attempts to minimize calculations by way of delayed calculations.
Whenever a control field is modified, DYNAMO marks that value as changed,
and the mathematical objects that depend on it (single-slice exponents, start-to-t and
t-to-end propagators etc.) as stale. Only when a specific stale object is needed
(such as propagators when computing the fidelity), are the calculations performed.
At this point DYNAMO attempts to perform the minimal number of matrix exponentiations and
multiplications to arrive at the desired result.

See cache.m for the nitty-gritty details.


License
=======

All the code is released under the Lesser GNU Public License (LGPL) 3.0,
except where explicitly stated otherwise. Everything else
(documentation, algorithms, etc.) is released under the Creative
Commons Attribution-ShareAlike 3.0 License. See "LICENSE.txt" for details.

In one sentence: Normal scientific community licensing philosophy ---
take this work, do whatever you like with it, publish any enhancements
you make back to the community, and don't forget to add a reference to
the original authors.


Authors
=======

Shai Machnes         2010-2012
Ville Bergholm       2011-2015
