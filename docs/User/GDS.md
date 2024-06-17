# GDS Exporter

The GDSII format is used in lithography. This is a better implementation of the Qiskit-Metal GDS renderer as:

- It eliminates potential creases or cracks forming due to floating-point errors
- Provides additional functionality such as boolean operations to add/combine layers

The basic usage is:

```python
from SQDMetal.Utilities.MakeGDS import MakeGDS

leGDS = MakeGDS(design)

#Optional boolean layer via a logical OR operation on layers 0 and 1
new_layer_ind = leGDS.add_boolean_layer(0, 1, "or")

leGDS.export('first.gds')
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

