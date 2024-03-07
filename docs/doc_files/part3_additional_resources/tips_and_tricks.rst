.. _tips:

Tips and Tricks around DART
============================

This section provides some tips and tricks in dealing with DART.

Visualizing molecules with concatenated xyz files
---------------------------------------------------

DART makes extensive use of concatenated xyz files, in which multiple molecules are stored in a single file. In the DART output, these files will always start with `concat_` and end on `.xyz`. Concatenated xyz files are a very handy output for the user because they make it very easy to browse through many molecules in order to doublecheck them. For example, the ``ligandfilters`` module saves concatenated xyz files for each filter and for the passed ligands, so that you can browse through the ligands and doublecheck visually if your filters work correctly.

Concatenated xyz files can be visualized quickly with the free software ``ase``, which is automatically downloaded with DART. If you are in the DART conda environment or you downloaded ase separately, you can use the following command in your terminal to visualize the concatenated xyz file:

.. code-block:: bash

    ase gui <path_to_xyz_file>

**Display bonds:** Per default, ``ase`` does not display bonds. To display bonds, you can press `Ctrl` + `B` in the ase viewer. Some versions of ase have white bonds which clash with the white background. In this case, you can change the default settings of ase such that it uses a black background. To do that, create a file under the path `~/.ase/gui.py` and add the following line:

.. code-block:: python

    gui_default_settings['gui_background_color'] = '#000000' # black background

 For more details, please refer to the ``ase`` documentation.
