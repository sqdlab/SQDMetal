.. _palace-sims:

Palace Simulations
==================

All the classes inherit:

- :class:`~SQDMetal.PALACE.Model.PALACE_Model_Base`
- :class:`~SQDMetal.PALACE.Model.PALACE_Model_Base_RF` (which inherits the above base class) in the case of RF simulations.

Thus, there refer to these classes for common interface functions.

To run capacitance simulations:

.. code-block:: python
   :emphasize-lines: 2

   from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

.. autoclass:: SQDMetal.PALACE.Capacitance_Simulation.PALACE_Capacitance_Simulation
   :members:
   :undoc-members:
   :show-inheritance:

To run RF eigenmode simulations:

.. code-block:: python
   :emphasize-lines: 2

   from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

.. autoclass:: SQDMetal.PALACE.Eigenmode_Simulation.PALACE_Eigenmode_Simulation
   :members:
   :undoc-members:
   :show-inheritance:

To run RF frequency driven simulations:

.. code-block:: python
   :emphasize-lines: 2

   from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

.. autoclass:: SQDMetal.PALACE.Frequency_Driven_Simulation.PALACE_Driven_Simulation
   :members:
   :undoc-members:
   :show-inheritance:
