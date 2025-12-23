.. _sqdmetal-api:

API Reference
=============

The SQDMetal API brings together different technologies 

Components from SQDMetal
------------------------

SQDMetal has custom components that are compatible with Qiskit-Metal


Simulations using AWS Palace
----------------------------

SQDMetal supports the ability to directly simulate models in Qiskit-Metal and GDS formats. The relevant constructs are:

- :ref:`Palace Simulation Classes <palace-sims>`
- :ref:`Palace Visualisation <palace-vis>`

.. toctree::
   :maxdepth: 2
   :hidden:

   palace
   palacevis


Useful Utilities
----------------

- :ref:`QUtilities <utilities-qutilities>` - a general set of helper utility functions involving Qiskit-Metal entities.
- :ref:`MakeGDS <utilities-gds>` - a GDS exporter with useful manipulation/postprocessing functionality.
- :ref:`ManipulateGDS <utilities-mgds>`
- :ref:`Qubit Designer <utilities-qbder>`

.. toctree::
   :maxdepth: 2
   :hidden:

   qutilities
   utilitygds
   utilitymgds
   qubitdesigner
