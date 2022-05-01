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

    # Find motifs in a single fasta file
    find_motifs_fasta("<fasta filepath>", "<output folder>")

    # Find motifs for multiple fasta files
    find_motifs_fasta_folder("<input folder>", "<output folder>")


# Find motifs for a single fasta file
"\<fasta filepath\>": a string that's the input fasta file's absolute filepath.

"\<output folder\>":  a string that's the output folder's absolute filepath; all the motif discovery results will be stored here.

# Find motifs for multiple fasta files
"\<input folder\>": a string that's the input folder's absolute file path; the input folder contains multiple fasta files.

"\<output folder\>": a string that's the output folder's absolute file path; this output folder will contain multiple folders. Each folder will store the motif discovery results that correspond to a input fasta file.


# Software requirements:

# Hardware requirements:

# Cite
More details on this later.

# Contact
Please contact 

