# Sphinx Documentation

## Building locally

The installation of SQDMetal already includes `sphinx`. Thus, enter the associated virtual environment (changing the name to the venv name) and run the following commands in the root directory of the cloned repository:

```bash
conda activate sqdmetal_env
sphinx-build doc-sphinx/source doc-sphinx/build
```

Then open `doc-sphinx/build/index.html`.

## Editing structure

It is recommended that one installs the extension `reStructuredText` by `lextudio` if using *VSCode*.

## Docstring convention

The convention was chosen to maintain readibility while being compatible with `sphinx`. Here is an example covering all the important details (note the important line breaks):

```python
"""
Comment on the function

Parameters
----------
param1 : str 
    Description of parameter un
param2 : int 
    Description of parameter deux
options : dict
    This is a dictionary (e.g. kwargs):

    *   ``'key1'`` (`float`):
            Description of this key
    *   ``'key2'`` (`int`):
            Description of this other key

Returns
-------
    Details of the return structure. Use similar structure to the dict above if returning dictionaries.
"""
```


