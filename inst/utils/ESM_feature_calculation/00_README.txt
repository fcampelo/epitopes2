Calculate ESM features

This series of functions/ routines is used to calculate
features using the ESM embedders, particularly ESM-1b.

Note that this uses Python3 and needs at least Pandas
and torch to be installed

Execute the files as follows:
- Edit 01_setup.R to configure your feature calculation parameters.
- source("01_setup.R")
- source("02_build_fasta.R")
- source("03_invoke_ESM.R")
 (or open the file and use Python commands in there to call directly from the terminal)
- source("04_concatenate_outputs.R")


IMPORTANT: This is a **very** computationally demanding activity,
and usually cannot run in regular computers (this function is intended
for use in high-performance workstations or clusters). The function also
needs a good internet connection to run, as the ESM model requires the
download of certain elements to be run.
To check the details and install instructions of ESM, please see
<https://github.com/facebookresearch/esm> and examples therein. A clone
of that repository is also available at <https://github.com/fcampelo/esm>.


ESM-1b returns errors when trying to calculate features for proteins longer
than 1024 residues (see https://github.com/facebookresearch/esm/issues/49).
A possible solution (suggested in, e.g.,
https://github.com/brianhie/evolocity/issues/2)
is to break longer sequences in smaller overlapping windows, predict
independently and then aggregate in post-processing, averaging where needed.
We follow this strategy here.

NOTE: For the feature calculation, any non-standard AA character is replaced by `-`, which is internally interpreted by ESM-1b as the `<unk>` token (see https://github.com/facebookresearch/esm/issues/13).

