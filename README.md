# Molecule Breakdown Model

Thank you for your interest in this repo complementing the "[XXX Analysis of the Generated Database GDB-13s](https:)" publication.

## Requirements
### You have installed anaconda or miniconda with python 3.6 - 3.9
### Install conda environment

The requirements for the environment are given in the fragment.yml file
<Br/>`conda env create -f fragment.yml`

Specific packages used are also listed below:
  - ipykernel >= 6.4.1
  - numpy >= 1.21.2
  - pandas >= 1.4.1
  - pickleshare >= 0.7.5
  - rdkit >= 2021.09.4
  - tqdm >= 4.63.0
  
### Conda environment activation
 `conda activate env-fragment`
 
## Quickstart

### Example

*example1.smi* and *example2.smi* files in the folder *example* are provided for the Molecule Breakdown Model demonstration.

  - The model can be executed by using the following command to obtain the ring fragments and corresponding substituents:

    `python Molecule_Breakdown_Model_Ring_Fragments.py example1.smi`
    
    `python Molecule_Breakdown_Model_Substituents.py example1.smi`

  - The model can also be carried out by applying the script *Molecule_Breakdown_Model_Ring_Fragments_Ring_Fragments.ipynb* and *Molecule_Breakdown_Model_Ring_Fragments_Substituents.ipynb* in the Jupyter Notebook

Two .pickle files containing the ring fragments and substituents for each molecule in the sample pool (you can find these data in folder *example/results_of_each_example*) will be obtained. The final ring fragment and substituent datasets after combining the duplicates and sorting their frequency (you can find these final fragment datasets for two examples in folder *results_merging/results* for further steps of merging, etc.). 

A .html file elaborating the final fragments results can be also obtained for each approach.

### Results Merging

For very large databases, we usually have to split the databases or even use high-performance computer clusters to treat them parallelly. 

Therefore, an efficient script for merging all the results for each sub-database is necessary. 

In folder *results_merging* you can see the file *Molecule_Breakdown_Results_Merge.ipynb* which will realize this merging process.

Then all sorted frameworks or the top 10000 sorted frameworks for the entire database can be displayed in a .html file. You can also save the SMILES of frameworks on your own or utilize the *results_merged.pickle* file, and then visualize them with software like ChemDraw or Marvin.

## Contributing

We welcome contributions, in the form of issues or pull requests.

If you have a question or want to report a bug, please submit an issue.

To contribute with code to the project, follow these steps:
1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Make your changes and commit them: `git commit -m '<commit_message>'`
4. Push to the remote branch: `git push`
5. Create the pull request.

### Contributors

* [@Ye-Buehler](https://github.com/Ye-Buehler)

The contributors have limited time for support questions, but please do not hesitate to submit an issue (see above).

## License

The software is licensed under the MIT license (see LICENSE file), and is free and provided as-is.
