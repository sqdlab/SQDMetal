# Overview of SQDMetal components

The components given in SQDMetal follow the conventions of *Qiskit-Metal* while adding some conventions to ensure uniformity. The names the components have hints as to their usage. For example, take `InductorMeander`:

- `InductorMeander` - positions itself via `(pos_x, pos_y)` and `(end_x, end_y)`
- `InductorMeanderPinStretch` - positions itself via a component pin and spans a distance `dist_extend` from this pin
- `InductorMeanderPinPin` - positions itself via across two component pins

Note the following:

- The components should import properly if SQDMetal was installed as instructed.
- Currently the *qlibrary* will not show the components in the *Qiskit-Metal* GUI toolbar. 
- However, the *Qiskit-Metal* GUI will show documentation of SQDMetal components (i.e. it will show the doc-string in its help section) and will have them behave correctly.
