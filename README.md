# CDLmotif

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/CDLmotif.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/CDLmotif.jl/dev)
[![Build Status](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/CDLmotif.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/CDLmotif.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/CDLmotif.jl)


# Introduction

This is a motif discovery method based on convolutional dictionary learning. More details will be posted later on.

# Usage:

    using CDLmotif

    # find motifs in a single fasta file
    # input: a string that's the input fasta file's absolute filepath
    # ouput: a string that's the output folder's absolute filepath; output folder will contain the motif discovery results on this fasta file
    find_motifs_fasta("<fasta filepath>", "<output folder path>")

    # multiple fasta files; find motifs in each fasta file
    # input: a string that's the input folder's absolute filepath; the input folder contains multiple fasta files)
    # ouput: a string that's the output folder's absolute filepath; output folders will contain a subfolder for each input fasta file in the input folder
    find_motifs_fasta_folder("<input folder that contains multiple fasta files>", "<output folder>")


# Software requirements:

# Hardware requirements:

# Cite
More details on this later.

# Contact
Please contact 

