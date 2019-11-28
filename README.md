# physics-toolbox

MRI modelling tools.

## Dependencies

Most functions implemented in this package depend on (a recent version 
of) the [SPM](https://www.fil.ion.ucl.ac.uk/spm/) software.

## Contents

### `+ismrmrd`

[ISMRMRD](https://ismrmrd.github.io) is a file format aimed to raw MR 
data (complex, k-space, multi-coil, non-processed). The `+ismrmrd` 
folder contains a reader for this format built around MATLAB's h5 
library.

### `+b1m`

The `+b1m` folder contains algorithms for the estimation of complex 
B1 receive fields from complex, multi-coil data. Its state is
work-in-progress.

### `+sense`

The `+sense` folder contains a very basic implementation of the SENSE 
algorithm for acceerated data acquired with a regular Cartesian sampling 
scheme. It assumes that complex coil sensitivities are known.

### `+mpm`

The `+mpm` folder contains algorithms for the estimation of multiple
parametric maps (R1, R2\*, PD, MT-sat) from multi-echo spoiled 
gradient-echo acquisitions. Its state is work-in-progress.

### `+optim`

The `+optim` folder contains generic numerical optimisation algorithms.

### `+utils`

The `+utils` folder contains various utility functions share by all 
above sub-projects.

## Contributors

This toolbox was developed at the [Wellcome Centre for Human 
Neuroimaging](http://www.fil.ion.ucl.ac.uk/), University College London.

If you'd like to try it out and encounter any difficulty, please send 
an email to `y.balbastre` *at* `ucl.ac.uk`

## License

This software is released under the [GNU General Public License version 
3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify 
the software as long as you track changes/dates in source files. Any 
modifications to or software including (via compiler) GPL-licensed code 
must also be made available under the GPL along with build & install 
instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))