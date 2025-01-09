# DART - Notes for developers

- The conda environment for development is `conda_envs/conda_DART_dev.yml`, see the README in that directory for more information.
- The documentation is built with Sphinx, see the README in the `docs` directory for more information.
- The tests are integration tests, which run an entire module and then check that the output is as expected.

## Installation

DART is installed via pip. One issue with DART is the installation of the MetaLig database, since this database is very big (393 MB). This file size is too big to upload for both GitHub and PyPI. Therefore, the MetaLig database is compressed via bzip2 and uploaded as compressed file. Once the installation is done, the database is decompressed and stored in the same directory as the compressed database.

**Potential issue:** The installation procedure was tested on Mac and seems to work, but it requires that python has the permission to write to the directory where the database is stored. This might be an issue on some systems, especially if people would install their python packages with sudo. However, right now I don't think it's a major issue, if it comes up we might have to install the database in a different way.

## Release

As preparation, make sure you installed twine and build, and added the PyPI and TestPyPI API credentials to your ``~/.pypirc`` file. These are helpful links:
* https://packaging.python.org/en/latest/guides/distributing-packages-using-setuptools/
* https://stackoverflow.com/questions/53122766/best-workflow-and-practices-for-releasing-a-new-python-package-version-on-github

Then follow these steps to release a new version of DART on PyPI (pip):
On your test branch, test on TestPyPI first:
   1. Append .1 to the end of the version number in ``setup.py`` and ``__init__.py`` (i.e. a 'debug' version ``D``) in case we need to upload multiple test versions while debugging.
   2. Build package locally: ``python3 -m build --sdist; python3 -m build --wheel``
   3. Upload to TestPyPI: ``twine upload -r testpypi dist/DARTassembler-X.Y.Z.D*``
   4. Make new conda environment to test the new version: ``conda create --name test_DARTassemblerX.Y.Z.D python=3.10 pip``. 
   5. Activate the new environment, then install and test from TestPyPI: ``pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple DARTassembler==X.Y.Z.D``
   6. If everything works, continue with the next steps. If not, fix the issue, increment the 'debug' version number and try again.
If everything works:
   7. Increment version number, this time properly:
      1. Set ``version=X.Y.Z`` in ``setup.py``
      2. Set ``__version__=X.Y.Z`` in package ``__init__.py``
   8. Push changes to master: ``git push`` with comment ``Bump to version X.Y.Z.``
   9. Build package locally again: ``python3 -m build --sdist; python3 -m build --wheel``
   10. Upload to PyPI Production: ``twine upload dist/DARTassembler-X.Y.Z*``
   11. On GitHub, create a new release with the tag ``vX.Y.Z`` and add a description of the changes.
   12. Re-build the documentation on ReadTheDocs: https://readthedocs.org/projects/dartassembler/builds/
   



