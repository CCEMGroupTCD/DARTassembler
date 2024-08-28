# DART - Notes for developers

- The conda environment for development is `conda_envs/conda_DART_dev.yml`, see the README in that directory for more information.
- The documentation is built with Sphinx, see the README in the `docs` directory for more information.
- The tests are integration tests, which run an entire module and then check that the output is as expected.

## Installation

DART is installed via pip. One issue with DART is the installation of the MetaLig database, since this database is very big (393 MB). This file size is too big to upload for both GitHub and PyPI. Therefore, the MetaLig database is compressed via bzip2 and uploaded as compressed file. Once the installation is done, the database is decompressed and stored in the same directory as the compressed database.

**Potential issue:** The installation procedure was tested on Mac and seems to work, but it requires that python has the permission to write to the directory where the database is stored. This might be an issue on some systems, especially if people would install their python packages with sudo. However, right now I don't think it's a major issue, if it comes up we might have to install the database in a different way.