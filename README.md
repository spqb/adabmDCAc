# adabmDCAc-

This is a C/C++ implementation of the Boltzmann machine learning for biological sequences.

## Installation
In the __src__ folder run
```
make
```
It will generate the executable file __adabmDCA__. See 
```
./adabmDCA -h
```
for a complete list of features.

## Basic runs


### Learning a Potts model from a MSA

```
./adabmDCA -f <MSA file> -k <label> -b <alphabet flag> -m <X iterations>
```
  - Output files will be saved every X iterations specified in the `-m` flag;
  - The output folder is named after the `label` specified in the `-k` flag.
  - Use `-w <file name>` for ad-hoc weights file (optional).

### Learning a sparse Potts model from a MSA

#### Decimation from converged run
```
./adabmDCA -f <MSA file> -k <label> -I <argument> -A -c <convergence threshold> -x <required sparsity>
```
  - `-A` flag removes gauge invariance at the beginning of the training;
  - `-I` allows the machine to start from a profile model for the fields, and the (rescaled by argument) covariance matrices for the couplings.

#### Activation from profile model
```
./adabmDCA -f <MSA file> -k <label> -I 0. -Z -c <convergence threshold> -x <required sparsity>
```
  - `-Z` flag inactivates all couplings at the beginning of the training.

### Restore training
Add the flag `--restore` to restart the training from the checkpoint saved in the output folder.

## Output files
  - A lot of things
  - _Parameters file_. Missing couplings are assumed to be inactive when the training re-starts.

  


