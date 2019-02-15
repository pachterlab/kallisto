---
layout: page
title: "Single-cell RNA-Seq analysis"
description: ""
---

The analysis of single-cell RNA-Seq data involves a series of steps that include: (1) pre-processing of reads to associate them with their cells of origin, (2) possible collapsing of reads according to unique molecular identifiers (UMIs), (3) generation of feature counts from the reads to generate a feature-cell matrix and (4) analysis of the matrix to compare and contrast cells.

Some of these challenges are procedurally straightforward but computationally demanding. Others are are statistical in nature and require technology specific models. We have recently introduced a format for single-cell RNA-seq data called the BUS (Barcode, UMI, Set) format that facilitates the development of modular workflows to address the complexities of these challenges. It is described in [P. Melsted, V. Ntranos and L. Pachter, "The Barcode, UMI, Set format and BUStools", bioRxiv 2018](https://www.biorxiv.org/content/early/2018/11/21/472571).

BUS files can be generated from single-cell RNA-seq data produced with any technology and can, in principle, be produced by any pseudoalignment software. We have implemented a command in __kallisto__ version [0.45.0](http://pachterlab.github.io/kallisto//releases/2018/11/17/v0.45.0) called "bus" that allows for the efficient generation of BUS format from any single-cell RNA-seq technology. Tools for manipulating BUS files are provided as part of the [__bustools__](https://bustools.github.io/) package. Finally, [R](https://github.com/BUStools/BUS_notebooks_R) and [python](https://github.com/BUStools/BUS_notebooks_python) notebooks for processing and analyzing BUS files simplify and facilitate the process of developing and optimizing analysis workflows.


{% include JB/setup %}


