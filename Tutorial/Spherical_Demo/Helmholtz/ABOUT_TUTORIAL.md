# Tutorial for Working on a Sphere {#tutorials3}
[TOC]

This tutorial walks through analysing a velocity field using Helmholtz decomposition with coarse-graining.

---

There are three stages in this tutorial, each with a corresponding sub-directory. See the ABOUT documentation file in each directory for the step-by-step guide to each part.
1. The 'project' directory works through the steps needed to perform the Helmholtz decompiosition
2. The 'filter' directory then illustrates how to produce spatial maps of the coarse-grained flow
3. The 'postprocess' directory shows how to compute area-averaged values and perform a scan of filter scales to produce a power spectrum

In this directory there is a python script `generate_data.py`.
This should be run before starting the other pieces.
The generated velocity field is a combination of randomly-seed vortices with specified size distribution.
 - full spherical domain, 0.5 degree resolution, polar 'continents'
 - velocity field is a collection of eddies of varying sizes, randomly placed throughout the globe

