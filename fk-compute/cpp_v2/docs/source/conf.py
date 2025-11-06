# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FK Computation C++ Library'
copyright = '2024, FK Computation Team'
author = 'FK Computation Team'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output ------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Todo settings
todo_include_todos = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'cpp': ('https://en.cppreference.com/w/', None),
}

# Math settings
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'

# Custom CSS
def setup(app):
    app.add_css_file('custom.css')

# Source file suffixes
source_suffix = {
    '.rst': None,
    '.md': 'myst_parser',
}

# Master document
master_doc = 'index'

# Language
language = 'en'

# Pygments style
pygments_style = 'sphinx'

# HTML theme options
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# HTML context
html_context = {
    'display_github': True,
    'github_user': 'fk-computation',
    'github_repo': 'fk-compute',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}

# HTML title
html_title = f'{project} v{release}'

# HTML short title
html_short_title = 'FK Computation'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_logo = None
html_favicon = None

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

# Additional templates that should be rendered to pages
html_additional_pages = {}

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer.
html_show_copyright = True

# If true, an OpenSearch description file will be output
html_use_opensearch = ''

# Output file base name for HTML help builder.
htmlhelp_basename = 'FKComputationdoc'