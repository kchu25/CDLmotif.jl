# CDLmotif

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/CDLmotif.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/CDLmotif.jl/dev)
[![Build Status](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/CDLmotif.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/CDLmotif.jl)


# Introduction

This is a DNA-sequence-motif discovery method based on convolutional dictionary learning. More details on the implementations and derivations of the model are coming soon. 

# Usage:
In Julia, import CDLmotif, and then use any of the two following subroutines:
```julia
using CDLmotif

# Find motifs in a single fasta file
# <fasta-filepath>: A string that's the input fasta file's absolute filepath.
# <output-folder>: A string that's the output folder's absolute filepath; 
#                  all the motif discovery results will be stored here.
find_motifs_fasta(<fasta-filepath>, <output-folder>)

# Find motifs for multiple fasta files
# <input-folder>: A string that's the input folder's absolute file path; 
#                 the input folder contains multiple fasta files. This folder 
#                 *must* contain only fasta files.
# <output-folder>: A string that's the output folder's absolute file path;
#                  this output folder will contain multiple folders. Each 
#                  folder will store the motif discovery results that 
#                  correspond to an input fasta file.
find_motifs_fasta_folder(<input-folder>, <output-folder>)
```

Please check below for a quick reference on the usage and software/hardware requirements.

# Find motifs for multiple fasta files

The input folder should look something like this:

    input-folder
    ├── <fasta_1>.fa
    ├── <fasta_2>.fa
    ├── <fasta_3>.fa
    ├── ...
    └── <fasta_K>.fa

The output folder will look something like this:

    output-folder
    ├── <fasta_1>
    |   ├── logos
    |   ├── ├── d1.transfac
    |   ├── ├── d2.transfac
    |   ├── ├── ...            
    |   └── summary.html
    ├── <fasta_2>
    ├── <fasta_3>
    ├── ...
    └── <fasta_K>

Motifs that are discovered will be stored in the *transfac* format as count matrices (for which it is straightforward to transform them into PWMs). Note that ```d1.transfac``` corresponds to the first discovered motif, and ```d2.transfac``` correspond to the second discovered motif, and so on. A summary on the motif discovery results is documented in ```summary.html``` in each folder.

# Data requirements:

### Fasta
We require the input data to be in [fasta format](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp). Currently, we only support fasta files that contain all sequence that have equal lengths (equal number of nucleotides in each sequence). 


# Software requirements:
### Weblogo
 This software currently requires the installation of [Weblogo](http://weblogo.threeplusone.com/manual.html#download). You can install Weblogo via the following command (with python3 installed beforehand):
 ```bash
 pip3 install weblogo
 ```
And you can check if you have weblogo installed by typing in your command line by the following command:
```bash
weblogo -h
```

### MEME
This software currently requires the installation of [MEME](https://meme-suite.org/meme/doc/download.html) (We use MEME's utility function fasta-shuffle-letters to create a control dataset to calculate the significance of each discovered motif using fisher exact test). Once MEME is installed, make sure you add MEME's utilities path to your PATH environment. On Linux operating systems, you can do so by adding the following line to your .bashrc:

```bash
export PATH=<path to where MEME is installed>/meme/bin:$PATH
```

# Hardware requirements:
 We require the user to have an Nvidia GPU as we currently use [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl/) to speed up the optimization for convolutional dictionary learning.

# Cite
More details on this later.

# Contact
Please contact <skchu@wustl.edu> or raise an issue on the github repo with any questions about installation or usage.

# Notes
- We plan to drop the dependence on MEME's utilities soon.
- Julia uses a just-in-time compiler, which means that software written in Julia needs to be pre-compiled before its execution. We hope to remove this soon so that we don't have to wait every time on pre-compilations before we actually run the subroutines. This makes the application more easily integrated into the bioinformatics pipeline (e.g., snakemake; see point below).
- We plan to extend this software to a standalone application so that the Julia language will no longer be required. This can be done with packages like [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl/). However, there's currently a bug related to [packages that depend on LLVM.jl](https://github.com/JuliaLang/PackageCompiler.jl/issues/682) that prevents this from being realized. We may implement a c++ version of this software if the LLVM-dependency issue is not resolved soon.
- Currently, we rank the motifs on the result page via the sum of likelihood ratio scores. We will add the motif significance for each discovered motif soon.
- More documentation and extensions on optional inputs on adjusting the hyperparameters are coming soon.

We have not added this package to the Julia registry yet. To use this software, one way to do so is to simply clone it, and enter the following in your julia code
```julia
push!(LOAD_PATH, <path to the folder that contains CDLmotif.jl>)    
```
before you type 
```
using CDLmotif
```
