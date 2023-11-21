.. _current_limitations:

Current Limitations of DART
=============================

In the latest version, DART still has some limitations which we will seek to address in future releases. If any of them are important for your work, please let us know and leave us a :ref:`feature request <contributing>`. These limitations include:

Assembler:
    - The hydride ligand (H) cannot yet be assembled due to issues with the underlying stk library.
    - Tridentate ligands with non-planar coordinating atoms cannot yet be assembled.
    - Only octahedral and square-planar geometries are supported. Let us know if you need other geometries for your research.

MetaLig database:
    - Ligands might have missing Hydrogen atoms or wrong charges. This should be extremely rare, but cannot be excluded because of the large number of ligands in the database.

Other limitations are more fundamental and will probably not be addressed anytime soon:

MetaLig:
    - Unique ligands are defined by their graph when connected to an (anonymous) metal center. Therefore stereoisomers are not distinguished in the ligand extraction process and only one stereoisomer would be part of the MetaLig database.

