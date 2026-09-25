"""viva-cpm — a process-bigraph Cellular Potts Model framework.

A fast Rust CPM engine (``viva_cpm.cpm_core``) with a process-bigraph Python layer:
CPM and subcellular processes, schema-driven world construction, connectivity /
basement-membrane structural constraints, and analysis metrics.

Typical use from another project::

    from viva_cpm import load_world, cpm_core            # engine + schema builder
    world = load_world(spec)                        # build a CPM world from a dict
    world.step(100)

Cells are wired as process-bigraph composites via import-path addresses, e.g.
``local:!viva_cpm.processes.cpm_process.CPMProcess`` and
``local:!viva_cpm.subcellular.sbml.SBMLSubcell`` — see ``cpm.composites`` for examples.
"""
from viva_cpm import cpm_core
from viva_cpm.schema import load_world

__all__ = ["cpm_core", "load_world"]
__version__ = "0.1.0"
