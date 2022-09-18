from .pangraph_ import PanGraph
#UMAP, JUMAPBASE

# import numba

import pkg_resources

try:
    __version__ = pkg_resources.get_distribution("PanGraph").version
except pkg_resources.DistributionNotFound:
    __version__ = "0.1-dev"