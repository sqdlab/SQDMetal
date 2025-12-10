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




