Qubit designer
==============

The final qubit and resonator system has a bunch of parameters following a set of design equations. The modules found under `SQDMetal.Utilities.QubitDesigner` facilitate in optimising said parameters. Currently, the supported resonator types are:

- λ/2 stripline resonator: :ref:`ResonatorHalfWave`
- λ/4 stripline resonator: :ref:`ResonatorQuarterWave`
- Lumped element LC resonator: :ref:`ResonatorLC`

The idea is that the qubit designer classes gobble up the resonator object to create the transmon-resonator system. Supported qubit types are:

- General transmon capacitively coupled to resonator: :ref:`XmonDesigner`
- Floating transmon coupled to resonator: :ref:`FloatingTransmonDesigner`

Resonator types
---------------

The resonator classes listed in this section all implement certain functions that can be useful:

- :func:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase.get_res_capacitance()` - returns the equivalent capacitance of the resonator (in Farad)
- :func:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase.get_res_inductance()` - returns the equivalent inductance of the resonator (in Henry)
- :func:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase.get_res_impedance()` - returns the equivalent impedance of the resonator (in Ω)
- :func:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase.get_res_frequency()` - returns the resonant frequency (in Hertz)
- :func:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase.is_res_parallelLC()` - returns whether the resonator is specified as a series or parallel LC circuit.


ResonatorHalfWave
~~~~~~~~~~~~~~~~~

This represents a λ/2 stripline resonator. The class calculates the equivalent capacitance and inductance depending on whether it is a open or short circuit termination:

.. code-block:: python

    from SQDMetal.Utilities.QubitDesigner import ResonatorHalfWave

    resonator = ResonatorHalfWave(7.25e9, shorted=False, impedance=50)

The first argument is the resonant frequency (in Hertz). The last two arguments are optional with the default values shown above; that is, a 50Ω λ/2 resonator with an open-ended termination (thus, equivalent to a parallel LC circuit).

ResonatorQuarterWave
~~~~~~~~~~~~~~~~~~~~

This represents a λ/4 stripline resonator. The class calculates the equivalent capacitance and inductance depending on whether it is a open or short circuit termination:

.. code-block:: python

    from SQDMetal.Utilities.QubitDesigner import ResonatorQuarterWave

    resonator = ResonatorQuarterWave(7.25e9, shorted=True, impedance=50)

The first argument is the resonant frequency (in Hertz). The last two arguments are optional with the default values shown above; that is, a 50Ω λ/4 resonator with a shorted termination (thus, equivalent to a parallel LC circuit).

ResonatorLC
~~~~~~~~~~~

This represents a simple parallel LC circuit. The class can be instantiated by providing any two of: inductance :math:`L`, capacitance :math:`C` and resonant-frequency :math:`f_0`. For example:

.. code-block:: python

    from SQDMetal.Utilities.QubitDesigner import ResonatorLC

    resonator1 = ResonatorLC(f0=7.25e9, L=1e-9)
    resonator2 = ResonatorLC(f0=7.25e9, C=1e-12)
    resonator3 = ResonatorLC(L=1e-9, C=1e-12)


Qubit types
-----------

The qubit designer classes all implement an ``optimise`` function:

- The function takes up a dictionary of parameters pertaining to that qubit type
- Each parameter in the dictionary can be either given as a single floating-point value (that is, a fixed constraint) or a 2-tuple (that is, a bounded interval).
- The ``optimise`` command will try to satisfy all the fixed and interval constraints. Note that if any of the parameters are not provided, a wide interval constraint is taken as a default value.

The function prints a list of the qubit parameters (e.g. :math:`g`, :math:`\chi`, anharmonicity etc.).

XmonDesigner
~~~~~~~~~~~~

This class is for a simple transmon capacitively coupled to a resonator:

.. code-block:: python

    from SQDMetal.Utilities.QubitDesigner import ResonatorHalfWave, XmonDesigner

    resonator = ResonatorHalfWave(10.5e9)
    XmonDesigner(resonator).optimise({
        'fQubit':9.5e9,
        'C_g':(0.1e-15,10e-15),
        'C_J':(10e-15,150e-15),
        'chi':(-10e6,-0.1e6),
        'Ej/Ec':(1,100)})


FloatingTransmonDesigner
~~~~~~~~~~~~~~~~~~~~~~~~

This class is for a floating transmon capacitively coupled to a resonator:

.. code-block:: python

    from SQDMetal.Utilities.QubitDesigner import ResonatorHalfWave, FloatingTransmonDesigner

    resonator = ResonatorHalfWave(10.5e9)
    FloatingTransmonDesigner(resonator).optimise({'fQubit':(1e9, 10e9), 
                                                'C_q1':(30.862e-15),
                                                'C_q2':(31.974e-15),
                                                'C_g1':(10.910e-15),
                                                'C_g2':(0.916e-15),
                                                'C_J':(42.461e-15, 55e-15),
                                                'chi':(-0.9e6,-0.2e6),
                                                'Ej/Ec':20})


