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

The above commands will require availability of atlas-orca data files. The location of atlas data files is by default:
 - For build-dirs: `<build-dir>/share`
 - For install-dirs: `<install-dir>/share`
 - Extra search paths can be given by `ATLAS_DATA_PATH` environment variable containing ':'-separated list of read-only search paths.

Data files can also be downloaded on-demand at runtime by setting following environment variables:
```bash
export ATLAS_ORCA_CACHING=1
export ATLAS_CACHE_PATH=<atlas-cache-dir>  # If not specified, `/tmp/cache` is used.
```
When a data file is not found in the search locations described above, they will be downloaded to the `ATLAS_CACHE_PATH`.


The data files can also be cached beforehand with a provided tool `atlas-orca-cache`. Use its `--help` argument for more info.

```bash
export ATLAS_ORCA_CACHING=1
export ATLAS_CACHE_PATH=<atlas-cache-dir>
atlas-orca-cache --grid=<gridname>  # use "all" as gridname to proceed for all available grids
```

You can then for following runs add this `ATLAS_CACHE_PATH` to the `ATLAS_DATA_PATH`.

Installing data files
---------------------

All orca data files can also be retrieved and installed as part of the build/install steps.
In order to proceed, configure the build with following CMake arguments:
```bash
-DENABLE_RETRIEVE_ORCA_DATA=ON -DENABLE_INSTALL_ORCA_DATA=ON
```
This could lead to several minutes of download times for each build. That is not a problem for central deployment.
However the use of a centralised ATLAS_DATA_PATH could be more beneficial, in a development environment with multiple users (e.g. prepIFS experiments).
