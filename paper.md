---
title: 'FlowSieve: A Coarse-Graining Utility for Geophysical Flows on the Sphere'
tags:
  - coarse graining
  - scale decomposition
  - C++
  - fluid dynamics
authors:
  - name: Benjamin A. Storer
    orcid: 0000-0001-5955-2158
    affiliation: 1
  - name: Hussein Aluie^[Corresponding author]
    orcid: 0000-0003-3516-3697
    affiliation: 1

affiliations:
 - name: University of Rochester, USA
   index: 1
date: 8 March 2022
bibliography: paper.bib

---

# Summary of Subject Area

Ocean and atmosphere dynamics span an incredibly wide range of spatial and temporal
scales, with spatial scales ranging from the sub-millimetre viscous scales all the way
up to planetary scales at tens of thousands of kilometres. Because of the strong non-linear
nature of oceanic and atmospheric flows, not only do the behaviour and characteristics change
significantly with scale, but important energetic interactions exist between different scales.
As a result, an important step in understanding and predicting the behaviour of such complex
flow systems is being able to disentangle the complex interactions across this wide range of scales.
Coarse-graining is a physically-motivated and mathematically-rigorous technique for partitioning
spatial flows as a function of a specified partitioning scale, allowing for a consistent and comprehensive
scale-by-scale analysis.

# Summary of Software

The core features of `FlowSieve` are:
1) computes coarse-grained scalar and vector fields for arbitary filter scales, in both Cartesian and spherical coordinates
2) built-in diagnostics for oceanographic settings, including kinetic energy (KE), KE cascades, vorticity, divergence
3) built-in post-processing tools compute temporal and region averages, for an arbitrary number of custom user-specified regions [ avoiding storage concerns when handling large datasets ]
4) includes Helmholtz-decomposition scripts to allow careful coarse-graining on the sphere [ i.e. to maintain commutativity with derivatives ]

`FlowSieve` is written in C++, with some python user-friendliness scripts included. 
Input and output files are netCDF.
`FlowSieve` is designed with heavy parallelization in mind, as well as several context-based optimizations, in order to facilitate processing high-resolution datasets. 
In particular, MPI is used to divide time and depth [ with minimal communication costs, since coarse-graining is applied at each time and depth independently ], while OpenMP is used to parallelize latitude and longitude loops, taking advantage of shared memory to reduce communication overhead.


`FlowSieve` can currently only work on mesh grids ( i.e. latitude grid is independent of longitude, and vice-versa ), but those grids need not be uniform.
This is not a restriction of the coarse-graining methodology, however, and future developments may extend this functionality if there is sufficient interest / need.


# State of the Field

Coarse-graining is being increasingly used as an analytical method in oceanographic communities. While the coarse-graining is similar to blurring / convolutions in image-processing, those tools do not readily apply to these contexts; they often rely on uniform, rectangular, Cartesian grids, which typically do not apply in Global Climate Model (GCM) data.
An established package in the fied, GCM-Filters, was designed to work on GCM data and grids.  It uses a diffusion-type coarse-graining method that is made available to users through python utilities.

The unique contributions of `FlowSieve` to the field are: follows the rigorous underlying mathematical framework to preserve physical properties of the data, designed for use on full spherical geometries, and can apply any arbitrary filtering scale spanning from sub-grid to domain-size.


# Statement of need

@Aluie2018 demonstrated how, when applied appropriately, coarse-graining can
not only be applied in a data-processing sense, but also to the governing equations.
This provides a physically meaningful and mathematically coherent way to quantify not
only how much energy is contained in different length scales, but also how much energy
is being transferred to different scales.

`FlowSieve` is a heavily-parallelized coarse-graining codebase that provides
tools for spatially filtering both scalar fields and vector fields in Cartesian
and spherical geometries. Specifically, filtering velocity vector fields on a sphere
provides a high-powered tool for scale-decomposing oceanic and atmospheric flows 
following the mathematical results in @Aluie2019.

`FlowSieve` is designed to work in high-performance computing (HPC) environments in order to
efficiently analyse large oceanic and atmospheric datasets, and extract scientifically meaningful
diagnostics, including scale-wise energy content and energy transfer.


# Acknowledgements

We would like to thank Mahmoud Sadek, Shikhar Rai, Michele Buzzicotti, Hemant Khatri, Stephen Griffies, for valuable feedback in developing FlowSieve.

Development of FlowSieve was financially supported by US NASA grant 80NSSC18K0772 and NSF grant OCE-2123496.

# References
