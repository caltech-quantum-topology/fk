# Configuration file for the Sphinx documentation builder.
#
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
    'sphinx.ext.mathjax',
    # 'breathe',
    # 'exhale',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Breathe Configuration ---------------------------------------------------

breathe_projects = {
    "fk-compute": "_build/doxygen/xml"
}
breathe_default_project = "fk-compute"

# -- Exhale Configuration ----------------------------------------------------

exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = ../include\n"
                             "RECURSIVE = YES\n"
                             "EXTRACT_ALL = YES\n"
                             "GENERATE_LATEX = NO\n"
                             "GENERATE_HTML = NO\n"
                             "GENERATE_XML = YES\n"
                             "XML_OUTPUT = _build/doxygen/xml\n"
                             "PREDEFINED = POLYNOMIAL_TYPE=0\n"
                             "ENABLE_PREPROCESSING = YES\n"
                             "MACRO_EXPANSION = YES\n"
                             "EXPAND_ONLY_PREDEF = NO\n"
                             "SEARCH_INCLUDES = YES\n"
                             "INCLUDE_PATH = ../include\n"
                             "FILE_PATTERNS = *.h *.hpp *.cpp\n"
                             "EXCLUDE_PATTERNS = */.git/* */build/*\n"
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

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

# Math support
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'

# Todo extension
todo_include_todos = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}