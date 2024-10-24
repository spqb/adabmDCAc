# adabmDCA 2.0 - Direct Coupling Analysis in C/C++

**Authors:**  
- **Lorenzo Rosset** (Sorbonne Université, Sapienza Università di Roma)
- **Roberto Netti** (Sorbonne Université)
- **Anna Paola Muntoni** (Politecnico di Torino)
- **Martin Weigt** (Sorbonne Université)
- **Francesco Zamponi** (Sapienza Università di Roma)
  
**Maintainer:** Anna Paola Muntoni

## Overview

**adabmDCA 2.0** is a flexible yet easy-to-use implementation of Direct Coupling Analysis (DCA) based on Boltzmann machine learning. This package provides tools for analyzing residue-residue contacts, predicting mutational effects, scoring sequence libraries, and generating artificial sequences, applicable to both protein and RNA families. The package is designed for flexibility and performance, supporting multiple programming languages (C++, Julia, Python) and architectures (single-core/multi-core CPUs and GPUs).  
This repository contains the C/C++ version of adabmDCA, maintained by **Anna Paola Muntoni**.

## Features

- **Direct Coupling Analysis (DCA)** based on Boltzmann machine learning.
- Support for **dense** and **sparse** generative DCA models.
- Available on multiple architectures: single-core and multi-core CPUs, GPUs.
- Ready-to-use for **residue-residue contact prediction**, **mutational-effect prediction**, and **sequence design**.
- Compatible with protein and RNA family analysis.

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
  - The output folder is named after the `output folder` specified in the `-a` flag. Files will be labeled according to the argument of the flag `-k`.
  - Use `-w <file name>` for ad-hoc weights file (optional).
  - For RNA, set the flag `-b n`; for ad hoc alphabet set `-b <alphabet>` where `alphabet` is a string

### Learning a sparse Potts model from a MSA

#### Decimation from converged run
```
./adabmDCA -f <MSA file> -k <label> -a <output folder> -p <params> -c <convergence threshold> -x <required sparsity> -L
```
  - `-A` flag removes gauge invariance at the beginning of the training;
  - Additional options `-V <drate>`

#### Activation from profile model
```
./adabmDCA -f <MSA file> -k <label> -a <output folder> -I 0. -Z -c <convergence threshold> -X <gsteps> -L -e <nsweep>
```
  - `-Z` flag inactivates all couplings at the beginning of the training;
  - `-I 0.` allows one to start from a profile model;
  - Additional options `-U <nactive>`;
  - Convergence at target Pearson whatever the density.

### Restore training
Add the flag `--restore` to restart the training from the checkpoint saved in the output folder.

### Sampling
Use
```
./adabmDCA -p <params> -f <MSA file> -i 0 -S -L -s <nconfig>
```
  - `-W nmix` (optional)


### Computing energies
```
./adabmDCA -p <params> -f <MSA file> --energies
```


## Output files
  - A lot of things
  - _Parameters file_. Missing couplings are assumed to be inactive when the training re-starts.

## To be done
  - Compute DMS score
