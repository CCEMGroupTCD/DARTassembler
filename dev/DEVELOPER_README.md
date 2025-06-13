# DART - Notes for developers

- The conda environment for development is `dev/conda_envs/conda_DART_dev.yml`, see the README in that directory for more information.
- The documentation is built with Sphinx, see the README in the `docs` directory for more information.
- The tests are integration tests, which run an entire module and then check that the output is identical to a benchmark directory.

## Refactoring for version 1.1 (rca_ligand_refactor branch):
- The `rca_ligand_refactor` branch is the `master` branch for the completely revamped version of DART, which will then be version 1.1.0.
    - Active developing should be done on other branches, which are then merged into the `rca_ligand_refactor` branch once the feature is complete and all tests pass.
      - Note: During the currently ongoing refactoring process, the installation, Pd-Ni and OER tests (see `Tests`) are not yet fully functional, because I believe it will make more sense to adapt them once at the end. All the other tests though are functional and should be used.


### Refactoring goals, plan and priorities
Note from Timo: In my time writing scientific python code, I have often tried to implement all the features I wanted from the beginning and to make everything perfect from the beginning on. However, every time I ended up with a very complex and slow development process, and in the end I still had to change a lot because scientific software development is always messy and one can never predict everything. On the other hand, whenever I restrained myself and focused on implementing just the core features first, I ended up with a much simpler and faster development process, and in the end I could implement all the advanced features by making only small, local changes. Therefore, I would like to learn from this experience now for the DART refactoring. From my experience I would say that it is very recommendable to develop in stages, where first one tries to develop the core functionalities asap and tests them well and completely forgets about the advanced stuff. Afterwards, the code will be in a good enough state that usually it is easy to go in and locally change and add something to implement advanced features, and it is much faster and better than when doing that from begin with. Therefore, I would like to suggest the following plan for the DART refactoring process:

Core functionalities to be implemented first, in order of priority:
1. Agree on the format for all the input of the `Isomer().from_ligands_and_metal_centers()` method which is used in the DART workflow. If any changes should be really required later on (which tbh will always be the case), let's talk and exchange our ideas early and agree on the updated specifications. For the addition of new features, new parameters can initially simply be added for quick development, but after a bit of a try-out phase we should exchange ideas and agree on the exact format of the new parameters, so we don't have to change it later on. 
2. Implement output of warnings for clashing ligands and symmetrically equivalent ligands (as can be seen at the end of the `Isomer().from_ligands_and_metal_centers()` method). Very important to do asap because it was part of the old DART workflow.
3. Implement everything to work perfectly with haptic and multi-metallic systems. This alone will take a good while and is difficult, so we need to focus on this because these are the actual core features that make DART so powerful. 
Lower priorities in the Isomer() class, to implement once the main priorities are fully done and well tested:
4. By default, remove the `try_all_geometrical_isomer_possibilities()` from `get_rotated_ligands()`, since as discussed, we will assume that the target vectors are in the correct order, and otherwise throw an error in the Python code for the user, or sort it ourselves in the DART workflow.
5. Add background atoms, an option to specify custom atoms which will be added to the complex xyz. This feature will probably also replace or be integrated with the old way of customized moving ligand atoms, as done in the rotation of the Pd-Ni ligand. 
6. Implement that the metal-ligand distance is adapted depending on the metal center a ligand is coordinated to. This should be done for now only for ligands which connect to only one metal center, not for bridging ligands, because that's much more complex. However, it is complicated, because it also needs to work for haptic ligands, where we have a Cu atom right now as the dummy atom, but of course the ligand should not be placed that far way, but instead depending on its actual donor atoms. There are two options I see how to implement this: keep track of the elements of the original haptic donor atoms and then calculate the mean donor distance from these, or changing the Cu dummy atom to one of the actual donor atom elements and then just use that as an approximation. I would personally prefer the first options because I think it is more robust and logical. 
7. Swap around ligands to be able to generate all isomers, also of monodentates. While this is a cool feature and I think it will definitely be part of DART soon, we are right now in the stage of developing and testing the core features, which is already difficult enough. This feature is quite complex and would slow down the development a lot if we take it into account all the time while developing the other stuff. Therefore, the development should happen in stages: In the first stage, we focus on the core stuff, which is already very difficult to get right. Then, once we have all the other code finished, it will be much easier to implement this feature as well.
8. Doublecheck and validate the user input. Very important but also time-consuming and prone to implement things that become unnecessary later on.
   - Small note about the last point, maybe it helps you: I have seen that you implemented very sophisticated classes for the input, like MetalSpec and LigandSpec. This is very professional a great idea.  However, right now we are in the stage of implementing the core features. Therefore, it is best to leave this doublechecking to the very end, once everything else is done, because it slows down the development of the actual features and usually just makes one write code that in the end isn't needed. I made exactly this mistake with the old DART workflow, where I spent a lot of time on this and then I had to change it all the time and now I can barely use any of it. Therefore, it's so important for me now to prioritize very hard, in order to learn from this experience. Additionally, I would also recommend to not use a yaml file but to simply define the input lists in the python file you use as playground for the development, since this is also easier and more intuitive (for an example see the end of the `assembler.py` file where I do this for the workflow). If you do prefer to have a yaml file (which in the end is equivalent), I would still recommend to not implement any doublechecking of the input for now.

## Code structure of the `DARTassembler` package
The DART workflow has several important Python classes that need to be maintained. One re-occurring feature is that for all DART modules that are accessible via the CLI, there should be a function to call this class from a yaml file. However, all classes themselves should take only Python objects (lists, dictionaries, strings, ase.Atoms etc) as input, since this makes the development much simpler. 

As such, the following Python classes need the ability to be accessed via the CLI and a yaml file:
- `Assembler()` -> the full DART workflow, from loading the ligand db files to saving the output 3D TMC structures
- `LigandFilters()` -> the ligand filters
- `DBInfo()` a module to get a csv and concatenated .xyz files from a ligand db .jsonlines file
- `Configs()` -> a module to generate a yaml file with the default configuration for the DART workflow, which can then be used as a template for the user to create their own yaml file
- `Concat()` -> a module to concatenate multiple ligand db files into one single file, which can then be used in the DART workflow
- BaseModule -> a base class for all modules that are accessible via the CLI, which implements the common functionality to read in a yaml file and call the module with the specified parameters. This is used by the Assembler, LigandFilters, Concat, Configs and DBInfo classes.

The following Python classes do not need to be accessed via yaml since they will never be called from the CLI
- `Isomer()` -> the class that represents 3D TMCs
- `Ligand()` -> the class representing a ligand from the MetaLig
- `IsomerFactory()` -> a class to generate Isomer() objects from a list of ligands and metal centers

### Assembler() 
This class runs the full DART assembly workflow and is maintained mostly by Timo. It implements the following functionalities:
- Run from a list of keywords (called in Python) or from a yaml file
- Load all required ligand databases
- Generate combinations of ligands to be assembled together in the same complex
- Call the `Isomer()` class to generate a list of 3D geometries from the specified ligands and metal centers
- Save all generated isomers and accompanying information to an output directory

### Isomer()
This class is used to assemble the 3D geometries of the ligands and is maintained mostly by Cian. It currently has the following inputs (and a few more options will still be added):
```python
class Isomer(BaseMolecule):
    def __init__(self,
                    atomic_props: Union[ase.Atoms, Dict[str, Any]],
                    graph: nx.Graph,
                    metal_idc: List[int],
                    donor_idc: List[List[int]],
                    ligand_idc: List[List[int]],
                    ligand_info: Dict[str, Any] = None,
                    global_props: Dict[str, Any] = None,
                    validity_check: bool = False,
                    ):
```
In the `Assembler()` class it is used like this:
```python
isomers, warnings = Isomer.from_ligands_and_metal_centers(
                                                                    ligands=ligands,
                                                                    target_vectors=self.target_vectors,
                                                                    ligand_origins=self.ligand_origins,
                                                                    metal_centers=self.metal_centers
                                                                    )
```
The input of this class are mostly Python lists of ligands, metal centers etc, which it uses to assemble a list of 3D geometries, which it returns together with a list of warnings (see above). It does not require to be run from a yaml file since it will not be called from the CLI, but instead it is a flexible and powerful Python class to be used for the general assembly of 3D transition metal complexes from a list of given ligands, metal centers, target vectors and origins (so it can also be used independently of the `Assembler()` class). It implements the following functionalities:
- Assemble a list of isomers from a set of ligands, metal centers, target vectors and origins
- Support for multi-metallic systems by associating ligands with each metal center in the input lists
- Support for haptic systems by treating a haptic group as a single effective donor atom situated in the center of the haptic group
- Check if two isomers are equivalent
- Check if ligands clash

### Ligand()
This class is used to represent a ligand from the MetaLig database and is maintained mostly by Timo. It has a wide range of properties to filter for in the ligandfilters module. It also has important properties for the `Isomer()` class such as the ligand geometry (`self.geometry`) and a list of list of indices of it's effective donor atoms, where the outer list are the possible (not-equivalent) possible orientations of the ligand (`self.geometric_isomers_hapdent_idc`).

### LigandFilters()
A class that loads a ligand database file and call be called from a yaml in order to filter down the ligand database based on certain filters.

### DBInfo()
A class that loads a ligand database file and can be called from a yaml in order to generate a csv file with the ligand information and a concatenated .xyz file with all ligands in the database.

## Tests

During the entire development process, the code should be continually tested to check if the output is as expected. The tests are integration tests, which run an entire module and then check that the output files are the same as in a benchmark directory. If they are not, the tests will print information about the differences. All tests are located in the `tests/integration_tests/` directory. The following tests are available:

Very fast (a few seconds):
- `test_cli.py`: A short tests of all modules using the cli. This test is very fast but touches every single module so it's very handy.
- `test_assembler.py`: Tests the `assembler` module which is the DART workflow, from reading in a yaml file, loading small sample ligand databases, assembling a few geometries and saving the output complexes. The tests here should be a little more expansive in the future, at the moment they are very fast as well.

Currently not functional during the refactoring:
- `test_Pd_Ni_example.py`: A test which runs the `ligandfilters` and `assembler` modules for the Pd-Ni case study.
- `test_OER_example.py`: A test which runs the `ligandfilters` and `assembler` modules for the OER database.

Note: The tests of ligandfilters, installtest and dbinfo that once existed have been removed because they are now all integrated in `test_cli.py`.

### How to test your code whenever you make any edits
When developing new features, it is important to test the code to ensure that the output is as expected. The tests are integration tests, which run an entire module and then check that the output files are the same as in a benchmark directory. If they are not, the tests will print information about the differences. This is a good way to ensure that the code is working as expected and that no unintended changes have been made.
As an example, let's imagine we want to develop a new feature in our assembler module. From experience, the following steps are a good practice to follow:
1. Make sure you have the latest version of the code by pulling (either from `master` or now from the `rca_ligand_refactor` branch).
2. Make sure you are in a git branch for testing and feature development.
3. Before making any edits, run all tests that contain the module we want to edit (or simply all of them) to make sure the benchmark outputs are correct. This will give us a baseline to compare against.
    - In our case, we edit the assembler module, so the most important test to run is `test_assembler.py`. 
    - The `test_cli.py` test should also always be run since it is very fast and covers all modules, including the assembler module.
    - The `test_Pd_Ni_example.py` and `test_OER_example.py` tests should also be run since the assembler module is used in these tests as well.
4. Check that all the tests pass successfully. This should normally be the case, but if not, first we need to fix the benchmark directories before we can continue. You should see the following output:
```
Integration test: check if the new output is the same as the old output.
Integration test successful: all good!
Test for ligand filters passed!
```
5. Now that we have a baseline, we can start editing the assembler module. Let's say we want to make a simple change such as changing the printed output, which should have no effect on the output files except the file `log.txt` which logs the output. So, knowing that, we will make the change and re-run the assembler test. The final output now looks like this:
```
Integration test: check if the new output is the same as the old output.
Changed: 1 item(s)
  /log.txt
==========    WARNING: INTEGRATION TEST FOUND ISSUES    ==========
	# changed files: 1
Test for assembly of complexes passed!
```
6. Now we need to check what has changed and if that is expected. The testing module will print a list of files which have changed. Here, we see that the `log.txt` file has changed, which is expected since we changed the printed output. We can also see that no other file has been changed because no other file is listed. So, depending on the complexity of the change we made, we can either continue now or we can use the `diff` functionality of PyCharm (see tips below) to further inspect the changes in the old and new `log.txt` file. Once we think that all changes are consistent with our expectations, we can continue with the next step.
7. If we are satisfied with the changes, we need to update the benchmark directory with the new output files in order to reflect the new gold standard output. For this, you simply go to the directory `tests/integration_tests/assembler/`. You find there three directories:
- `data_input`: This directory contains the input files for the test, which should not be changed.
- `data_output`: This directory contains the output files of the test, which have just been generated. 
- `benchmark_data_output`: This directory contains the benchmark output files, which were used to compare against the new output files.
Now, you can simply delete the `benchmark_data_output` directory and rename the `data_output` directory to `benchmark_data_output`. That way, the next time you run the test, it will make a new `data_output` directory and compare against the new `benchmark_data_output`.
8. Run the test again to make sure everything is working as expected. You should see the following output, indicating that the test passed successfully without any changes:
```
Integration test: check if the new output is the same as the old output.
Integration test successful: all good!
Test for ligand filters passed!
```

### Best practice: interactively testing your code during development
Above, I have described how to test the code and make sure that the output is as expected. However, this should usually not only be done once at the end, but rather continuously during the development process. This is a good practice to ensure that the code is working as expected at every step of the way. It helps to catch bugs early and makes it easier to debug the code since you can always see exactly the impact the code changes have on the output files. Apart from the Pd-Ni and OER example tests, the tests are made to be very fast, so you can run them interactively while developing the code. 

This approach is especially powerful if you make changes in the code that should not change the output files, such as refactoring, since you can immediately spot any issues. However, also when developing new features, I personally use the tests continually. Usually, I try to divide the development process into small steps, each of which can be tested individually. This way, I can run the tests after each step and see if the output is as expected. If not, I can immediately debug the code and fix the issue. Only once I'm satisfied with the output, I update the benchmark directory if necessary and continue with the next step. This way, I have a very strong control over the changes I introduce. 

Initially, this testing process might seem like a lot, but one gets used to it very quickly, and it will become a very powerful feeling to have such tight control over the output files and to always be able to spot any issue right from the bat. 

### Further tips for testing and debugging
- The `diff` functionality of PyCharm is very powerful and I found it way to late. In the project tab, you can mark two files and then right-click to select "Compare Files". This will show you the differences between the two files in a very nice way, highlighting the changes.


## Installation

DART is installed via pip (`pip install DARTassembler`). One issue with DART is the installation of the MetaLig database, since this database is very big (360MB). Therefore, the pip package includes a compressed version of the database (`MetaLigDB_v{VERSION}.jsonlines.bz2`) which is only 38MB big, which fits with the PyPI size limit of 50MB. The code that reads in the file will automatically detect if it's compressed and will uncompress it on-the-fly, line-by-line. Note that the database is never decompressed and written to disk because writing to the directory where something is installed after the installation is not a good practice and can lead to issues such as permission errors. Instead, the database is read directly from the compressed file. Only in the first version of DAT (up until version 1.1.0), the database was uncompressed and written to disk once at the beginning, which is not the case anymore.

## Release of a new DART version on PyPI

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

## Changelog
12.06.2025
- Removed the total charge and metal oxidation state from assembler.yml input and replaced this with a single property `total_ligand_charges`
- Refactored the tests. There is now a new test test_cli.py which tests all modules via the CLI. The old tests dbinfo, ligandfilters and installtest have been removed since they are now all integrated in the test_cli.py test. It is recommended from now on to use the cli test to always test everything, and if required also use the assembler test. The assembler test is right now very small, but it should be expanded in the future to test more features of the assembler module.
- Added a deprecation warning if test_metalig is used. From now on, every function and module has a property called either n or n_max_ligands, which is the maximum number of ligands to be read in, which should be used instead for testing purposes.
- Refactored the CLI and made a new general class BaseModule. This class is used by the Assembler, LigandFilters, Concat, Configs and DBInfo classes to provide a common interface for the CLI.
- Renamed the old ligand MetaLig ligand database files in data/metalig by prepending OLD_ to the file name to avoid confusion with the new ligand database files. The old files are still there for reference, but they should not be used anymore.


## Todo list
### Small renaming and shifting
- Consistently name the `n_max_ligands`. Right now often it is `n` in all modules except the assembler.
- rename complex_name_appendix to complex_name_suffix in assembler.yml
### Functionalities
- Add that in assembler.yml, the user can specify or leave out total_ligand_charges=None, which will then mean that any ligand combination is fine, independent of their total charge.