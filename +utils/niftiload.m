function vol = niftiload(fname)
% Loads a matlab array from a Nifti file
%
% FORMAT vol = utils.niftiload(fname)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

nii = nifti(fname);
vol = nii.dat();