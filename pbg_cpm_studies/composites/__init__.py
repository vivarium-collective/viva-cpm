"""Composite factories for cpm-studies.

Importing this package registers every ``@composite_generator``-decorated
composite (so the dashboard's discovery + the study-ref availability lint see
them). Keep new composite modules imported here.

Convention — composite vs. single-process config: a composite *composes* (wires
two or more processes, or adds emitters/structure). A "composite" that wraps a
single process with a hardcoded config is really a config preset of ONE process,
not a composition. For a family of configs of the same process, register ONE
parameterized composite and pass params from each study, e.g. ``recruitment`` in
``chemotaxis.py`` (parameters ``cue_rate`` / ``chemo_lambda``) — not one composite
per condition. ``scripts/lint-workspace.py`` warns on single-process composites
that carry no parameters. (The gg1993 composites predate this and are the reason
the lint currently warns; they could be collapsed into one parameterized composite.)
"""
from . import gg1993  # noqa: F401  (fires @composite_generator decorators)
from . import chemotaxis  # noqa: F401
