import os
import sys
sys.path.insert(0, os.path.abspath('../../SQDMetal'))  # Path to your project

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'nbsphinx'
]

# Automatically generate summary pages (tables) for classes, functions, and methods
autosummary_generate = True  # This will create individual .rst files for each class and function
add_module_names = False    #It improves the readibility of the documentation by truncating module name from documented items.
autoclass_content = 'both'
toc_object_entries_show_parents = 'hide'    #Stop it showing class-name prefix on navbar (c.f. https://github.com/pradyunsg/furo/discussions/692)

# Additional autodoc settings (optional)
autodoc_default_options = {
    'members': True,           # Include all members (functions, methods, etc.)
    'undoc-members': True,     # Include members without docstrings
    'private-members': False,  # Exclude private members (those starting with '_')
    'show-inheritance': True,  # Show inheritance for classes
}

# Autosummary settings
# autosummary_imported_members = True  # Include members of imported modules

# -- HTML output ------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'  # You can change this to any theme you like
html_static_path = ['_static']
html_theme_options = {
    'navigation_depth': 10
}

nbsphinx_execute = 'never'