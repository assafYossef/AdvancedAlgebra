# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath('..'))
project = 'Galwa Fields'
copyright = '2024, Ido Nahum and Assaf Yossef'
author = 'Ido Nahum and Assaf Yossef'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
names_to_skip = ['inverse_matrix', '_adjugate_matrix', '_transpose', '_cofactor_matrix', 'determinant', 'pad_element',
                 'valid_repr', 'zero_element_check', 'same_field', 'same_prime_field', 'refactor_polynom_terms']


def maybe_skip_member(app, what, name, obj, skip, options):
    if name in names_to_skip:
        return True
    return None


def setup(app):
    app.connect('autodoc-skip-member', maybe_skip_member)
