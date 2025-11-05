# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # Point to DARTassembler root directory

project = 'DARTassembler'
copyright = '2024, CCEM group'
author = 'Timo Sommer, Cian Clarke, Felix Kleuker, Max García-Melchor'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import moldoc.molecule as molecule

extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.autodoc',  # <-- reads docstrings from Python modules
    'sphinx.ext.autosummary',  # <-- auto-generates stub pages for each class/function
    # 'sphinx.ext.napoleon',  # <-- parses Google/NumPy-style docstrings
    'sphinx.ext.viewcode',  # <-- adds "View Source" links
    'sphinx_rtd_theme',
    'sphinx_toolbox.confval',
    'sphinx_design',
    'sphinxarg.ext',  # for documenting command-line interfaces
    "sphinx_copybutton",
    'moldoc',  # for rendering 3D interactive molecular structures in the documentation
    # 'nbsphinx',  # for including Jupyter Notebooks in the documentation
]

# Autodoc settings
autodoc_typehints = "none"   # show type hints in the parameter list
autodoc_member_order = 'bysource'   # preserves order of methods as in source code
autosummary_generate = True         # enables automatic stub generation
autoclass_content = 'both'          # 'both', 'class' or 'init'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# nbsphinx settings
nbsphinx_allow_errors = False  # fail if notebook execution errors
nbsphinx_execute = 'auto'      # execute only if outdated

# Set the default molecule configuration for moldoc
moldoc_default_molecule_config = molecule.MoleculeConfig(
    background_color=molecule.Color(252, 252, 252),    # match exact background color of the docs
    material=molecule.MeshPhongMaterial(),                              # set the material of the molecule, this looks best
    is_outlined=False,                                                  # outline looks bad
)

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

from docutils import nodes
from docutils.parsers.rst import roles

# Define custom role 'filter' to format filter names in code style
def filter_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    node = nodes.literal(text, text, classes=["filter"])
    return [node], []

def arg_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    node = nodes.literal(text, text, classes=["filter"])
    return [node], []

roles.register_local_role('filter', filter_role)
roles.register_local_role('arg', arg_role)

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/DART_pic.png'


# User defined
import os
import sys
sys.path.insert(0, os.path.abspath("../"))
# Add the version and release information to the documentation
version = release = os.environ.get('READTHEDOCS_VERSION', 'dev')