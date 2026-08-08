"""Composite factories for cpm-studies.

Importing this package registers every ``@composite_generator``-decorated
composite (so the dashboard's discovery + the study-ref availability lint see
them). Keep new composite modules imported here.
"""
from . import gg1993  # noqa: F401  (fires @composite_generator decorators)
from . import chemotaxis  # noqa: F401
