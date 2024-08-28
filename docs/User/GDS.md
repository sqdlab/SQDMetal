# GDS Exporter

The GDSII format is used in lithography. This is a better implementation of the Qiskit-Metal GDS renderer as:

- It eliminates potential creases or cracks forming due to floating-point errors
- Provides additional functionality such as boolean operations to add/combine layers
- Adds export options for certain lithographic processes (i.e. positive or negative)
- Allows export of certain layers
- Add rendered text as a seperate layer on the exported GDS

The basic usage is:

```python
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
```

Note that the returned layer index can be used for further boolean operations. The geometry can also be rounded in post-processing:

```python
from SQDMetal.Utilities.MakeGDS import MakeGDS

#Set the initial components to be with no rounding - can set it to have some rounding if desired (in units of metres)
leGDS = MakeGDS(design, curve_resolution=60, precision=1e-12, smooth_radius=0)    #Use 60pts per quarter rotation and ensure enough precision to have nice rounded corners

#Set the rounding radius to be non-zero now... Now all subsequent layers will have a rounding of 200nm
leGDS.smooth_radius = 200e-9
bottom_layer = gdsexp.add_boolean_layer(0,1,'or')
top_layer = gdsexp.add_boolean_layer(2,2,'or')

leGDS.export('rounded.gds')
```

### Adding text to the export

You can also add text with the `add_text()` function as follows:

```python
from SQDMetal.Utilities.MakeGDS import MakeGDS

leGDS = MakeGDS(design)

#Add a text label
leGDS.add_text(text_label="Add your label here", layer=10, size=600, position=(0.1, 0.9))

#Export design with text label on layer 10
leGDS.export('design_with_text.gds')
```
