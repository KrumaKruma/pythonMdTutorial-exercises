# Python MD Tutorial

In this tutorial we will learn Python in the context of computational physics. We will do two exercises which do both the same, a Molecular Dynamics (MD) simulation of Lennard-Jones clusters. In the first exercise you will do the simulation only using the numpy package (https://numpy.org/). In the second exercise we will use the more advances atomic simulation environment to do the MD simulation (https://wiki.fysik.dtu.dk/ase/). Please, solve all TODOs in the python code step by step as described below. Before you start some installations have to be done. 

## Installation
In case you are working with Wondows system you first have to install **WSL** via the windows store. This will give you a Linus environment on Windows. Once you have installed WSL you can start the installation the same way as on Linux.
First we install the (mini-)anaconda:
```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
```
and accept auto start conda envrionment. In a next step we check if conda is installed in the right version:
```bash
    conda --version
```
Once conda is installed we can create a conda environment:
```bash
    conda create --name PythonTutorial python=3.10
```
Next we activate the conda environment we just created so that we work in the right environment:
```bash
   conda activate PythonTutorial 
```
Now we can install the packages nessecairy to solve the exercises i.e. numpy and ase:
```bash
    conda install numpy
    conda install -c conda-forge ase
```
Next we look if the packages are installed correctly:
```bash
    python -c "import numpy; print(numpy.__version__)"
    python -c "import ase; print(ase.__version__)"
```
Now that we have installed all the nessecairy libraries we can start with the exercises.

## Exercise 1: MD implementation with numpy
Here we will write an MD using the velocity Verlet algorithm only using the numpy library. To be successfull you have to implement the following things:
1. Write a read function to read an xyz file (numpyMD/readWrite.py)
2. Create a structure object (numpyMD/md.py)
3. Write an initialize velocity function (numpyMD/boltzmannVelocity.py)
4. Initialize velocities according to temperature (numpyMD/md.py)
5. Implement the velocity Verlet algorithm (md.py)

## Exercise 2: MD implementation with ase
In this exercise we will implement an MD using the velocity Verlet algorithm using the atomic simulation environment. To solve this task successfully you should implement the following things:
1. Read the xyz-files using ase io functionality
2. Crate the Lennard-Jones calculator and attach it to the atoms object
3. Initialize velocities according to the temperature using ase
4. Implement the velocity Verlet algorithm using ASE

## Check your implementation
In order to check if your implementation was successfull you can look if the energy is conserved. Plot the potential and the total energy in order to verify that your implementation was successfull. You can also download the ovito-visualizer (https://www.ovito.org/) to look at the movements of the atoms during the MD.
