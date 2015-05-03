---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

Typing `kallisto` produces a list of usage options, which are:

    Usage: kallisto <CMD> [arguments] ..

    Where <CMD> can be one of:

        index         Builds a kallisto index 
        quant         Runs the quantification algorithm 
        h5dump        Converts HDF5-formatted results to plaintext
        version       Prints version information

    Running kallisto <CMD> without arguments prints usage information for <CMD>

The four usage commands are:

#### index

`kallisto index` builds an index from a Fasta formatted file of target sequences. 

#### quant

`kallisto quant` runs the quantification algorithm. 

#### h5dump

`kallisto h5dump` converts [HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to plaintext. 

#### version

`kallisto version` displays the current version of the software.