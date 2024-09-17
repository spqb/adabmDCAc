# adabmDCAc

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
./adabmDCA -f <MSA file> -a <output folder> -k <label> -m <nsave> -L
```
  - Output files will be saved every `nsave` iterations specified in the `-m` flag;
  - The output folder is named after the `output folder` specified in the `-k` flag.
  - Use `-w <file name>` for ad-hoc weights file (optional).
  - For RNA, set the flag `-b n`; for ad hoc alphabet set `-b <alphabet>` where `alphabet` is a string

### Learning a sparse Potts model from a MSA

#### Decimation from converged run
```
./adabmDCA -f <MSA file> -k <label> -p <params> -c <convergence threshold> -x <required sparsity> -L
```
  - `-A` flag removes gauge invariance at the beginning of the training;
  - Additional options `-X <drate>`

#### Activation from profile model
```
./adabmDCA -f <MSA file> -k <label> -a <output folder> -I 0. -Z -c <convergence threshold> -x <required sparsity> -L -e <nsweep>
```
  - `-Z` flag inactivates all couplings at the beginning of the training;
  - `-I 0.` allows one to start from a profile model;
  - Additional options `-U <nactive>`, `-X <gsteps>`.

### Restore training
Add the flag `--restore` to restart the training from the checkpoint saved in the output folder.

## Output files
  - A lot of things
  - _Parameters file_. Missing couplings are assumed to be inactive when the training re-starts.

## To be done
  - Compute DMS score, energies of a MSA
  - Output/Input including alphabet in parameter file
  - Sampling at convergence


