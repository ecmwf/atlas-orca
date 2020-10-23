atlas-orca
==========

This is a Atlas plugin that adds support for ORCA grids and mesh generation


Requirements:
-------------
- atlas 0.22.0 or greater
- git lfs

NOTE! Before cloning this repository, make sure git is recent enough to support git lfs,
      and you have executed at least once for your user account

    git lfs install


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

    atlas-meshgen ORCA1_T --generator="orca" --3d -o mesh.msh
