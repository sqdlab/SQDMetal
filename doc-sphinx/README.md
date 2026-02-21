# Sphinx Documentation

## Building locally

The installation of SQDMetal already includes `sphinx`. Thus, enter the associated virtual environment (changing the name to the venv name) and run the following commands in the root directory of the cloned repository:

```bash
conda activate sqdmetal_env
sphinx-build doc-sphinx/source doc-sphinx/build
```

Then open `doc-sphinx/build/index.html`.

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

Note the line-spacing required before the first outlined dictionary argument.

## Editing structure

It is recommended that one installs the extension `reStructuredText` by `lextudio` if using *VSCode*. To write in *reStructured Text* (RST), see the primer [here](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html) and the tips below:

### Sectioning

Use the underlining to create sections, but don't nest too much:

```rst
Main Title
==========

Some text

Subsection
----------

More text

Subsection
~~~~~~~~~~

Even more text

```

### Linking within document

Because the `sphinx.ext.autosectionlabel` is enabled, the link can be done via the title alone:

```rst
See :ref:`A Plain Title`.

With :ref:`custom text <Another title>`
```

### Linking to other documents

Use the path (can be relative or from source):

```rst
Look :doc:`here <simpalacedv>`
```

### Links to functions

To link to other functions:

```rst
Here is my link to a function :func:`SQDMetal.PALACE.Model.PALACE_Model_Base.prepare_simulation`
```

To only show the function name, use a `~`:

```rst
Here is my link to a function :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.prepare_simulation`
```

To have custom text:

```rst
Here is my link to a function :func:`custom text <SQDMetal.PALACE.Model.PALACE_Model_Base.prepare_simulation>`
```

### Links to websites

To link to a website:

```
Here is a link `Link text <https://www.example.com/>`__
```

Note that it must have a space before the `<` and there must be two trailing underscores `__`


### Code blocks

Code blocks look like (**note the mandatory line-space**):

```rst
.. code-block:: python

    np.array([3,1])
```

while inline code requires double ticks:

```rst
Here is some ``parameter`` written down.
```

### LaTeX

Use math blocks:

```rst
:math:`\lambda/2*e^x`
```

