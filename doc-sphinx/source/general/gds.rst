GDS operations
==============

This documentation covers the basic functionality of GDS-related operations in ``SQDMetal``, namely :ref:`exporting <GDS Exporter>`, :ref:`adding text to a design <Adding text to the export>`, and :ref:`manipulating an existing GDS file <GDS Manipulator>`.


GDS Exporter
------------

The GDSII format is used in lithography. :class:`SQDMetal.Utilities.MakeGDS` is a better implementation of the Qiskit-Metal GDS renderer as:

- It eliminates potential creases or cracks forming due to floating-point errors
- Provides additional functionality such as boolean operations to add/combine layers
- Adds export options for certain lithographic processes (i.e. positive or negative)
- Allows export of certain layers
- Add rendered text as a seperate layer on the exported GDS

The basic usage is:


.. code-block:: python

    from SQDMetal.Utilities.MakeGDS import MakeGDS

    leGDS = MakeGDS(design)

    #Optional boolean layer via a logical OR operation on layers 0 and 1
    new_layer_ind = leGDS.add_boolean_layer(0, 1, "or")

    #Export all layers
    leGDS.export('first.gds', export_type="all")

    #Flatten all layers and export the positive pattern
    leGDS.export('first_positive.gds', export_type="positive")

    #Export only layer 0
    leGDS.export('first_layer0.gds', export_layers=[0])

Note that the returned layer index can be used for further boolean operations. The geometry can also be rounded in post-processing:

.. code-block:: python

    from SQDMetal.Utilities.MakeGDS import MakeGDS

    #Set the initial components to be with no rounding - can set it to have some rounding if desired (in units of metres)
    leGDS = MakeGDS(design, curve_resolution=60, precision=1e-12, smooth_radius=0)    #Use 60pts per quarter rotation and ensure enough precision to have nice rounded corners

    #Set the rounding radius to be non-zero now... Now all subsequent layers will have a rounding of 200nm
    leGDS.smooth_radius = 200e-9
    bottom_layer = leGDS.add_boolean_layer(0,1,'or',output_layer=20)
    top_layer = leGDS.add_boolean_layer(2,2,'or', output_layer=21)

    #Export with boolean layers on 20 and 21
    leGDS.export('rounded.gds')


Adding text to the export
-------------------------

You can also add text with the :func:`~SQDMetal.Utilities.MakeGDS.MakeGDS.add_text()` function as follows:


.. code-block:: python

    from SQDMetal.Utilities.MakeGDS import MakeGDS

    leGDS = MakeGDS(design)

    #Add a text label in the top left of the design
    leGDS.add_text(text_label="Add your label here", layer=10, size=600, position=(0.1, 0.9))

    #Export design with text label on layer 10
    leGDS.export('design_with_text.gds')


GDS Manipulator
---------------

You can import and array an existing GDS design using the :class:`SQDMetal.Utilities.ManipulateGDS` class as shown in the following example (be sure to replace the relevant path arguments -- ``import_path_gds`` and ``export_path_gds`` -- with your own paths):

.. code-block:: python

    from SQDMetal.Utilities.ManipulateGDS import ManipulateGDS

    #Import the existing GDS `import_design.gds`, choose the export file, and shift the origin of the imported file.
    f = ManipulateGDS(
                    import_path_gds="/path_to_design/import_design.gds",
                    export_path_gds="/path_to_design/export_design.gds",
                    origin=(0, 4000)
                    )

    # make a 16 x 16 array of the design from `import_design.gds` with 500 um spacing in the array. Here we only import the cell 'TOP' from `import_design.gds`. The array will export a single flattened cell 'Array' to the path given in the ManipulateGDS() initialisation above.
    m = f.make_array_onChip(columns=16,
                            rows=16,
                            export=True, 
                            spacing=("500um", "500um"),
                            use_cells=['TOP']
                            )

    #TODO: Remove text below once ManipulateGDS has been documented.

The arguments for initialisation of a ``ManipulateGDS()`` object are: 

- ``import_path_gds`` - Path to GDS file for import
- ``export_path_gds`` - Path to export GDS
- ``import_cells`` - (Optional) Cells to import (must be cells present in the passed GDS), defaults to all
- ``origin`` - (Defaults to (0,0)) Tuple containing relative origin to import the design


The arguments for ``make_array_onChip()`` are:

- ``columns`` - number of columns in the array
- ``rows`` - number of rows in the array
- ``spacing`` - (Defaults to ("50um", "50um")) Tuple containing x- and y-spacing of the array as strings with units
- ``chip_dimension`` - (Defaults to ("20mm", "20mm")) Dimension of the chip on which to array the structure (WIP)
- ``export`` - (Defaults to True) Choose whether to export the generated file (default), or just to operate on the design
- ``export_path`` - (Optional) Export path (if different from self)
- ``use_cells`` - (Optional) List of cells from the input GDS to include in the arrayed structure
