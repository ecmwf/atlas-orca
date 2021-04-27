atlas-orca
==========

This is a Atlas plugin that adds support for ORCA grids and mesh generation


Requirements:
-------------
- atlas 0.24.0 or greater

Optional:
- NetCDF, in order to generate new ORCA grid data files (see tool atlas-orca-convert)

Use:
----
Add to environment:

    export ATLAS_PLUGINS=atlas-orca

If the plugin is installed in the same location as atlas, it will be automatically discovered
upon loading atlas during `atlas::initialize(...)`.

Otherwise, make sure the root dir is found at runtime by setting:

    export ATLAS_PLUGINS_SEARCH_PATHS=/path/to/atlas-orca_ROOT

Verify the plugin works:

    atlas-grids ORCA1_T --info

To generate a atlas-orca mesh:

    atlas-meshgen ORCA1_T --3d -o mesh.msh

