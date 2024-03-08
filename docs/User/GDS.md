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

Note that the returned layer index can be used for further boolean operations.
