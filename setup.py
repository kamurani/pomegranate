import io
import os
from setuptools import setup, find_packages

VERSION = None
HERE = os.path.abspath(os.path.dirname(__file__))
NAME = "POMEGRANATE"

DESCRIPTION = "Phosphosite motif explorer -- graph network abstraction through embeddings"

# Import the PYPI README and use it as the long-description.
# Note: this will only work if "README.md" is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))

setup(
    name="POMEGRANATE",
    version='0.1.0',
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'graphein',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'pomegranate = pomegranate.cli:main',
            'pom = pomegranate.cli:main',
            'po = pomegranate.cli:main',
        ]
    }
)
