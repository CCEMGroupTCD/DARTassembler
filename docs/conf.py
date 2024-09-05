# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DART'
copyright = '2024, CCEM group'
author = 'Timo Sommer, Cian Clarke, Felix Kleuker'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import moldoc.molecule as molecule

extensions = [
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
    'sphinx_toolbox.confval',
    "sphinx_copybutton",
    'moldoc',  # for rendering 3D interactive molecular structures in the documentation
]

# Set the default molecule configuration for moldoc
moldoc_default_molecule_config = molecule.MoleculeConfig(
    background_color=molecule.Color(252, 252, 252),    # match exact background color of the docs
    material=molecule.MeshPhongMaterial(),                              # set the material of the molecule, this looks best
    is_outlined=False,                                                  # outline looks bad
)

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/DART_pic.png'


# User defined
import os
import sys
sys.path.insert(0, os.path.abspath("../"))
