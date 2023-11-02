# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = u"danrerlib"
copyright = u"2023, Ashley Schwartz"
author = u"Ashley Schwartz"

# -- General configuration ---------------------------------------------------
def skip_private(app, what, name, obj, skip, options):
    if name.startswith('danrerlib.mapping._'):
        skip == True
    return skip

def skip_member(app, what, name, obj, skip, options):
    if "danrerlib." in name:
        if obj.name.startswith('_'):
            skip = True
        if obj.name.endswith('_path') or obj.name.endswith('_dir'):
            skip = True
        if obj.name.endswith('_URL') or obj.name.startswith('GO_PATH'):
            skip = True
        if obj.name.endswith('_URL') or obj.name.startswith('GO_IDS_PATH'):
            skip = True
    if name in ["danrerlib.utils", "danrerlib.settings", "danrerlib.database"]:
        skip = True
    return skip

def skip_util_classes(app, what, name, obj, skip, options):
    if "util" in name:
       skip = True
    return skip

def skip_subpackages(app, what, name, obj, skip, options):
    if what == "package":
        skip = True
    return skip

def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_member)
    # sphinx.connect("autoapi-skip-member", skip_subpackages)
    sphinx.connect("autoapi-skip-member", skip_util_classes)


extensions = [
    "myst_nb",
    "autoapi.extension",
    # "sphinx.ext.napoleon",
    # "sphinx.ext.viewcode",
]
autoapi_dirs = ["../src/danrerlib"]
autoapi_template_dir = "_templates/autoapi"
stylesheet = "_static/css/parsing.css"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages. 
html_theme = "renku"
html_static_path = ['_static']
html_css_files = ["css/custom.css"]
html_show_sphinx = False
html_logo = '_static/img/danrerlib_logo.png'
html_title = "danrerlib"
