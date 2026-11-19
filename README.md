# `transittraj` Paper Supplementary Code

## Introduction

Today's public transit vehicles produce a large amount of automatic
vehicle location (AVL) data. This data is very useful for planning and
performance studies, but can be noisy, error-prone, and sparse. 
We present `transittraj`, an R package which cleans AVL point data and reconstructs
continuous, differentiable, monotonic, and invertible vehicle trajectory functions.

<img src="vignettes/figures/figure_1.png" alt="Example `transittraj` trajectory." width="80%" />

This repository provides code to support the `transittraj` paper
(*under review*), including all data cleaning, analysis, and visualization code.
Unfortunately, we do not have permission to share the underlying data, though
we encourage the reader to replicate our analyses with their own data.
The vignettes presented here are not intended as stand-alone articles;
for a beginner-friendly introduction to `transittraj` with open sample data,
check out the vignettes on the
[package website](https://utel-uiuc.github.io/transittraj/).

## Navigating this Repository

We've split the code used in the paper
into four `markdown` vignettes, all stored in [`/vignettes`](vignettes/):

  1. [The AVL Cleaning Workflow](vignettes/cleaning.md). This vignette includes
  the code for the visualizations presented in Section 3
  (Package Architecture and Workflow) and the cleaning results presented in
  Section 4 (Case Study).
  
  2. [Estimating Signal Performance Measures Using Reconstructed Trajectories](vignettes/tspm.md).
  This vignette includes all code used to estimate and visualize the
  traffic signal performance measures presented in Section 4 (Case Study). Also
  included is a sensitivity analysis of threshold values, shown in Appendix C.
  
  3. [Evaluating Data Quality](vignettes/evaluation.md). This vignette
  includes all code used to the quality of AVL data and
  trajectories before and after the cleaning workflow, as discussed in Section
  5 (Data Quality Evaluation).
  
  4. [Comparison of Monotonic Enforcement Methods](vignettes/mono.md). This
  vignette includes all code used to compare the data quality of forward-
  stepping and the proposed bidirectional approaches, presented in Appendix A.
  
Each `markdown` file contains the code used, the output (figures, tables, and text)
of that code, and annotations describing the code. Thorough discussions of the
methods and results are presented in our paper. In addition to the rendered ("knit") files, each
vignette has a un-rendered `Rmarkdown` file in [`/vignettes`](vignettes/).
These include all raw code chunks, including those hidden in the knit
vignettes (loading data, saving results, rendering plots, etc.).
