"""
MolSSI_FSHXF
A Python package to study laser driven non adiabatic dynamics with mixed quantum-classical methods based on the Exact Factorization approach
"""

# Add imports here
from .functions import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
