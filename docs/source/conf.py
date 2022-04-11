# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))
# print(os.path.abspath('../..'))

#import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'skeleton'
copyright = '2022, jcschindler01'
author = 'jcschindler01'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  'sphinx.ext.autodoc',
  'sphinx.ext.mathjax',
  # 'sphinx_rtd_theme',
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# HTML Theme Options for sphinx_rtd_theme
html_theme_options = {
    #'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    #'analytics_anonymize_ip': False,
    #'logo_only': False,
    'display_version': False,
    'prev_next_buttons_location': 'none',
    #'style_external_links': False,
    #'vcs_pageview_mode': '',
    'style_nav_header_background': '#3ea8a8',
    # Toc options
    #'collapse_navigation': True,
    #'sticky_navigation': True,
    #'navigation_depth': 4,
    #'includehidden': True,
    'titles_only': True,
}

# -- Other Settings ----------------------------------------------------------

# Logo
html_logo = '_static/logo.png'


# # -- AutoDoc Settings ------------------------------------------------------
autoclass_content = 'both'
autodoc_member_order = 'bysource'
autodoc_docstring_signature = True
autodoc_typehints = 'both'


# # -- Old 'classic' theme settings -------------------------------------------
# html_theme = 'classic'
# Change sidebar. Method from stackexchange (https://stackoverflow.com/questions/18969093/how-to-include-the-toctree-in-the-sidebar-of-each-page).
# html_sidebars = { '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html'] }
