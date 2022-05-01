# CDLmotif

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/CDLmotif.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/CDLmotif.jl/dev)
[![Build Status](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/CDLmotif.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/CDLmotif.jl)


# Introduction

This is a motif discovery method based on convolutional dictionary learning. More details will be posted later on.

# Usage:

import the CDLmotif library

    using CDLmotif

find motifs in a single fasta file

input: a string that's the input fasta file's absolute filepath

ouput: a string that's the output folder's absolute filepath; output folder

    find_motifs_fasta("<fasta filepath>", "<output folder path>")

find motifs for multiple fasta files

input: a string that's the input folder's absolute file path; the input folder contains multiple fasta files)
 
output: a string that's the output folder's absolute file path; output folders will contain a subfolder for each input fasta file in the input folder

    find_motifs_fasta_folder("<input folder>", "<output folder>")

# Software requirements:

# Hardware requirements:

# Cite
More details on this later.

# Contact
Please contact 

