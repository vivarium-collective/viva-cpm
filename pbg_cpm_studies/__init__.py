"""pbg_cpm_studies — workspace Python package.

Importing the package registers the workspace's composites so the dashboard's
composite discovery and the study-ref availability check can see them.
"""
from . import composites  # noqa: F401  (registers @composite_generator composites)
