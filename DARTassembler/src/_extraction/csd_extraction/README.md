#CSD Extraction Scripts

This is the collection of all the scripts I have used to extract the data from the CSD. The repo is structured quite easily

## data

This is the folder to dump all the data in


## src

Herein you find all the python files.

I provided one example, how to generally extract a specific entry from the CSD, if you want to play around with the
objects you can extract from the CSD.

Moreover, I have some utility files, which are mostly not used anymore. They were mainly conversion between different
objects of different packages, I kept them for the sake of completeness.

Lastly, there are files starting with "Script_...", these are the actual extraction scripts and pretty straight forward to
read.


## Get started

1. Make sure, that the CSD software is installed and that you have activated it using the key

2. Add the following line to your .zshrc
    ```
    export CSDHOME=/Applications/CCDC/CSD_2022
    ```
    If you didn't select the standard installation path, you may have to change that path
    This is required for the API to work

3. Create the environment
    ```
    conda env create -f environment.yml
    ```

4. Now, you're ready to go.





