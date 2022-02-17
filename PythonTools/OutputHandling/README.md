/page OutputHandling1
# Output Handling

---

Coarse-graining a large dataset over a range of filter scales can produce a cumbersome amount of data, since it retains full space-time dependence and just adds a new filter-scale dimension.

Because of this, on-line postprocessing was introduced to the code that allows the user to output processed results (time averages, (custom) area averages, zonal averages, and others), which requires substantially less storage space.

That is, the on-line postprocess allows for a much more extensive set of filter scales, provided that full space-time dependence of the coarse-grained fields is not needed.

In order to allow the user to run to coarse-graining routine multiple times one different sets of filter scales, a different output file is created for each filter scale.

While this makes it easy to add new filter scales to the results without having to re-run to entire process, it can make it tedius to compare the results across filter scales.

## Result-stitching Scripts

To files are included here that make stitching the separate postprocessing results into one file with an filter-scale dimension.

### merge_postprocess_results.py

This is a python routine that takes in some command-line arguments that specify filename information ( see 'python merge_postprocess_results.py --help' ) and then stitches the postprocess files into one since netcdf file.

The filter-scale dimension (called ell) is prepended in all cases, so that it is always the first dimension. The script automatically sorts by filter-scale when merging, so there is no need to worry about the ordering of the filenames that are passed in.



### merge_postprocess_files.sh

This bash file simply runs the corresponding python script and illustrates sample usage.


## Downsample-stitching Scripts

Another useful technique is to downsample the original dataset when coarse-grainin at very large filter scales.
However, this once again produces results in different files, and the previous postprocess-stitching scripts will not work, since there are different grids involved.

A second set of routines is included, when combines the results from different down-samplings.

Simply first run the postprocess-stitching routine on the output files for each grid, and then pass the resulting files to the downsample-merging scripts.

### merge_downsampling_results.py

This is a python routine that takes in some command-line arguments that specify filename information ( see 'python merge_downsampling_results.py --help' ) and then stitches the different resolution files into one since netcdf file.

The script automatically determines which grid has the highest resolution.

**Note**: Duplicated filter scales are *okay*, with only the results from the highest resolution grid being retained.

**Note**: Where necessary, downsampled results are interpolated onto the highest resolution grid, using a nearest-neighbours method to preserve the 'coarseness'.

### merge_resolutions.py

This bash file simply runs the corresponding python script and illustrates sample usage.
