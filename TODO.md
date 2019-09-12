# TO DO

## Functionality

- [] Add option to do slice-by-slice loading. That is, only load in one time-point at a time (per mpi process). This would increase IO, but would greatly reduce the memory requirements in the event of a large dataset. This would be a little bit of a headache, and I suspect there will be some surprise issues lurking in the corners (specifically in regards to the postprocessing suites). My guess is that it would realistically take two or three days to get all of the bugs worked out and have a clean way of interfacing with the on-line postprocessing.


### Postprocessing
- [] Right now the post-processing produces a separate file for each filter scale. This is somewhat cumbersome and not entirely useful. If they were all in one file, with a filter-scale dimension, then it would be easier to do cross-scale comparisons. It would also save memory in the event that we write the integration regions to the postprocessing file. Conceptually, this should be straight-forward enough, it would just be a little tedius, since we would need to handle the case where the filter scales aren't sorted. Not hard, just needs some careful coding. This would probably take the better part of a day to go from start to finished / tested.
- [] Output the time-averaged fields. The main hold-up here is that, when Ndepth > 1, we may have also split depth across processors. In that case, time averaging isn't just a simpe reduction. We would need to compare depth regions to lineup the processor appropriately before reducing the arrays. This itself wouldn't be difficult, it would just take some care. The cleanest way would probably be to create a sub-communicator, which would be the group of processors that share a specific depth-slice. That way would could simply call the normal reduction operators, but it would only apply to the sub-comm object, not the full communication world. Again, this wouldn't be a difficult task, it would just take a bit of time to make sure it's done properly. Probably a day to get this done, start to finished / tested

## Unit Tests
- [] Add tests for differentiation tools
 - [x] spherical deriv on spherical grid
 - [] Cartesian deriv on spherical grid
 - [x] Cartesian deriv on Cartesian grid
- [x] tests for vorticity tools
