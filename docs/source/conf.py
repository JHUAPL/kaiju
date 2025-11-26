# Configuration file for the Sphinx documentation builder.
import os
# -- Project information

project = 'kaiju'
copyright = '2025, JHU/APL and NSF NCAR'
author = 'Kaiju Development Team'

release = '0.75.3'
version = '0.75.3'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

html_static_path = ['_static']
html_logo = '_static/MAGE_Logo_final_dark-bg_vertical.png'

html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'collapse_navigation': False,
    'navigation_depth': 4,
}

html_css_files = [
    'css/sidebar_theme.css',
]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

#if not os.environ.get("BUILD_ALL"):
#    exclude_patterns.append("_obsolete/**")

def setup(app):
    if os.environ.get('BUILD_ALL'):
        app.tags.add('obsolete')