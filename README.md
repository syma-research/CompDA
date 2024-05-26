This R package implements the CompDA method (compsitional differential abundance) for testing _conditional_ microbial differential
abundance with health outcome, as described in the paper  "Compositional Differential Abundance Testing: Defining and Finding a New 
Type of Health-Microbiome Associations" (Ma, Huttenhower, and Janson, 2024). 

* The package can be installed in R with `devtools::install_github("syma-research/CompDA")`.

* The main funciton is `CompDA`. You can read more on how to use in R it with `?CompDA`.

* A companion paper repository is available at (https://github.com/syma-research/CompDA_paper).

* This package is a protoype. We plan to expand on its user interface and functionality in the near future. We will also provide a
  vignette to showcase CompDA's application in practice. For now, a helpful example where you can see how CompDA is used can be
  found in the [data analysis script of the paper repo](https://github.com/syma-research/CompDA_paper/blob/main/mds/5.0_CRC.qmd),
  where we apply CompDA to study gut microbes that are conditionally differentially abundant in colorectal cancers.

* Parallelization is implemented in CompDA with the `future` framework. By default, CompDA will perform single-core computing.
  Parallel computation can be conducted by specifications with `future` options such as `future::plan(future::multisession())`.
  You can read more on how to specify parallelization using `future` on their [documentation page](https://www.futureverse.org/packages-overview.html).
