# transittraj_paper-supp

## Introduction

Today's public transit vehicles produce a large amount of automatic
vehicle location (AVL) data. This data is very useful for planning and
performance studies, but can be noisy, error-prone, and sparse. 
We present `transittraj`, an R package which cleans AVL point data and reconstructs
continuous, differentiable, monotonic, and vertible vehicle trajectory functions.

<img src="vignettes/figures/figure_1.png" alt="Example `transittraj` trajectory." width="80%" />

This repository provides code to support the `transittraj` paper
(*under review*), including all data cleaning, analysis, and visualization code.
Unforatunately, we do not have permission to share the underlying data, though
we encourage the reader to replicate our analyses with their own data.
The vignettes presented here are not intended as stand-alone articles;
for a beginner-friendly introduction to `transittraj` with open sample data,
check out the vignettes out our
[package website](https://utel-uiuc.github.io/transittraj/).

## Navigating this Repository

We've split the code used in the paper
into three `markdown` vignettes, all stored in [`/vignettes`](vignettes/):

  1. [The AVL Cleaning Workflow](vignettes/cleaning.md). This vignette includes
  the code for the visualizations presented in Section 2
  (Package Architecture and Workflow) and the cleaning results presented in
  Section 3 (Case Study).
  
  2. [Estimating Signal Performance Measures Using Reconstructed Trajectories](vignettes/tspm.md).
  This vignette includes all code used to estimate and visualize the
  traffic signal performance measures presented in Section 3 (Case Study).
  
  3. [Evaluating the Robustness of Reconstructed Trajectories](vignettes/robustness.md).
  This vignette includes all code used to quantify and visualize the error in
  reconstructed trajectories across various polling frequencies and interpolating
  methodologies, as discussed in Section 4 (Robustness Check).
  
Each `markdown` file contains the code used, the output (figures, tables, and text)
of that code, and annotations describing the code. Thorough discussions of the
results are presneted in our paper. In addition to the rendered files, each
vignette has a un-rendered `Rmarkdown` file in [`/vignettes`](vignettes/).
These include all raw code chunks (including those hidden in the final
vignettes). Finally, all high-resolution raster (PNG) versions of the final
graphics are saved in [`/vignettes/figures`](vignettes/figures/).


