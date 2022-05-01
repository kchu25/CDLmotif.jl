# CDLmotif

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/CDLmotif.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/CDLmotif.jl/dev)
[![Build Status](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/CDLmotif.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/CDLmotif.jl)


# Introduction

This is a motif discovery method based on convolutional dictionary learning. 

# Usage:
```julia
using CDLmotif

# Find motifs in a single fasta file
find_motifs_fasta(<fasta-filepath>, <output-folder>)

# Find motifs for multiple fasta files
find_motifs_fasta_folder(<input-folder>, <output-folder>)
```

Please refer below for a quick reference on the usage and software/hardware requirements.

# Find motifs for a single fasta file
```julia
# <fasta-filepath>: A string that's the input fasta file's absolute filepath.
# <output-folder>: A string that's the output folder's absolute filepath; 
#                  all the motif discovery results will be stored here.
```


# Find motifs for multiple fasta files
```julia
# <input-folder>: A string that's the input folder's absolute file path; 
#                 the input folder contains multiple fasta files. The folder 
#                 *must* contain only fasta files.
```
The input folder will look like this:

    input-folder
    ├── fasta_1.fa
    ├── fasta_2.fa
    ├── fasta_3.fa
    ├── ...
    └── fasta_K.fa

```julia    
# <output-folder>: A string that's the output folder's absolute file path;
#                  this output folder will contain multiple folders. Each 
#                  folder will store the motif discovery results that 
#                  correspond to an input fasta file.
```
The output folder will look like this:

    output-folder
    ├── fasta_1
    |   ├── logos
    |   └──  summary.html
    ├── fasta_2
    ├── fasta_3
    ├── ...
    └── fasta_K

# Software requirements:
This software currently requires installation of [MEME](https://meme-suite.org/meme/doc/download.html) (We use fasta-shuffle-letters to create a control dataset to calculate the significance of each discovered motif using fisher exact test). Once MEME is installed, make sure you add MEME's utilities path to your PATH environment. On Linux operating systems, this is done by adding the following line to your .bashrc:

```bash
export PATH=<path to where MEME is installed>/meme/bin:$PATH
```

# Hardware requirements:
 We require the user to have an Nvidia GPU as we currently use [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl/) to speed up the optimization for convolutional dictionary learning.

# Cite
More details on this later.

# Contact
Please contact <skchu@wustl.edu> or raise an issue on the github repor with any questions about installation or usage.

# Notes
- We plan to drop the dependence to MEME's utilitiies in the near future.
- Julia uses a just-in-time compiler, which means that software written in Julia needs to be pre-compiled before its execution. We hope to get rid of this so that this application can be more easily integrated into the bioinformatics pipeline (e.g., snakemake; see point below)
- We plan to extend this software to a standalone application so that the Julia language will no longer be required. This can be done from packages like [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl/). However, there's currently a bug related to [packages that depend on LLVM.jl](https://github.com/JuliaLang/PackageCompiler.jl/issues/682) that prevents this from being realized. We may implement a c++ version of this software if the LLVM-dependency issue doesn't get resolved soon.
