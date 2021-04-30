atlas-orca
==========

This is a Atlas plugin that adds support for ORCA grids and mesh generation.
The plugin registers known ORCA grids and is responsible for retrieving and 
reading binary data files that store the coordinates for these grids.



Requirements:
-------------
- atlas 0.25.0 or greater



Loading of atlas-orca plugin:
-----------------------------

- If your project explicitly links with the atlas-orca library, then nothing needs to be done.
  The plugin will be loaded as any explicitly linked library.

- If atlas-orca is installed in the same install-prefix as the executable or the eckit library, 
  then nothing needs to be done. The plugin will be automatically detected and loaded at runtime.
  This is the recommended approach.
  
- Otherwise, two environment variables need to be set, which will help detection and loading
  of plugins:

      export PLUGINS_MANIFEST_PATH=/path/to/atlas-orca/share/plugins  # colon separated list
      export LD_LIBRARY_PATH=/path/to/atlas-orca/lib:$LD_LIBRARY_PATH # use DYLD_LIBRARY_PATH for macOS
    
Verify the plugin is loaded using the `atlas` executable:

    atlas --info
    
should print a "Plugin" section containing "atlas-orca" with associated version and git sha1.


Using atlas-orca plugin:
------------------------

Standard tools that link with atlas only should now be able to use Atlas concepts that have been
implemented in the atlas-orca plugin. 

Examples are Grid registration and MeshGenerator

Verify the plugin works:

    atlas-grids ORCA1_T --info

To generate a atlas-orca mesh:

    atlas-meshgen ORCA1_T --3d -o mesh.msh

The above commands will download required data files containing grid coordinates etc. 
at runtime as needed into the Atlas "cache" directory (default=`/tmp/cache`). This directory can be controlled with
the environment variable: 

    export ATLAS_CACHE_PATH=/path/to/cache

The files can also be cached beforehand with a provided tool to prevent unexpected downloads

    atlas-orca-cache [--help]
    
Downloads at runtime can also be explicitly prevented with

    export ATLAS_ORCA_DOWNLOAD=0
