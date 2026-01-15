# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import glob
import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

pyproject_path = os.path.abspath("../pyproject.toml")
pyproject_data = {}

if os.path.exists(pyproject_path):
    with open(pyproject_path, "rb") as f:
        pyproject_data = tomllib.load(f)

poetry_tool = pyproject_data.get("tool", {}).get("poetry", {})


project = "ScanpEx"
copyright = "2026, Yuji Okano"
author = "Yuji Okano"
version = poetry_tool.get("version", "0.1.0")
release = version

dependencies = poetry_tool.get("dependencies", {}).keys()

exclude = {"python"}

package_map = {
    "scikit-learn": "sklearn",
    "biopython": "Bio",
    "opencv-python": "cv2",
    "pillow": "PIL",
    "pyyaml": "yaml",
    "scikit-image": "skimage",
}

mocks = []
for pkg in dependencies:
    if pkg in exclude:
        continue
    import_name = package_map.get(pkg, pkg)
    mocks.append(import_name)

autodoc_mock_imports = mocks

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    # 'nbsphinx',
    # 'sphinx_gallery.load_style',
    "myst_parser",
]

# Napoleon settings
napoleon_google_docstring = False
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


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
htmlhelp_basename = "scanpex_doc"
html_logo = "_static/scanpex_logo_large.png"
html_theme_options = {
    "navigation_depth": 5,
    "logo_only": True,
    # "sidebarbgcolor": "#003F67"
}
master_doc = "index"
latex_documents = [
    (master_doc, "scanpex.tex", "ScanpEx Documentation", "Yuji Okano", "manual"),
]

man_pages = [(master_doc, "scanpex", "ScanpEx Documentation", [author], 1)]

texinfo_documents = [
    (
        master_doc,
        "scanpex",
        "ScanpEx Documentation",
        author,
        "scanpex",
        "ScanPy Extension and kwarg Preferences",
        "Miscellaneous",
    ),
]

# nbsphinx_thumbnails = {
#     "/".join(
#         v.split(".")[:-1]
#     ): v.replace(
#         "notebooks", "_static"
#     ).replace(
#         "ipynb", "png"
#     ) if os.path.exists(
#         v.replace(
#             "notebooks", "_static"
#         ).replace(
#             "ipynb", "png"
#         )
#     ) else "_static/scanpex_logo_mini.png" for v in glob.glob("notebooks/*.ipynb")
# }

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

autodoc_typehints = "description"
