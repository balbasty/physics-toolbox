function out = extrapolate(out,suffix,pth0)
% Extrapolate log-linear fit to TE=0.
% Equivalently, take the exponential of the log-intercepts.
%
% FORMAT out = mpm.estatics.extrapolate(out,[suffix],[path])
% out    - Output structure
% suffix - Filename suffix [.nii]
% path   - Output path [next to input]

if nargin < 3
    pth0 = '';
end
if nargin < 2
    suffix = '.nii';
end 

if suffix(1) ~= '.' && suffix(1) ~= '_'
    suffix = ['_' suffix];
end

if isfield(out, 'PDw')
    out.PD0 = out.PDw;
    if isa(out.PD0.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.PD0.dat.fname);
        end
        out.PD0.dat.fname = fullfile(pth, ['PD0' suffix]);
        create(out.PD0);
    end
    out.PD0.dat(:,:,:) = exp(out.PDw.dat());
end
if isfield(out, 'T1w')
    out.T10 = out.T1w;
    if isa(out.T10.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.T10.dat.fname);
        end
        out.T10.dat.fname = fullfile(pth, ['T10' suffix]);
        create(out.T10);
    end
    out.T10.dat(:,:,:) = exp(out.T1w.dat());
end
if isfield(out, 'MTw')
    out.MT0 = out.MTw;
    if isa(out.MT0.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.MT0.dat.fname);
        end
        out.MT0.dat.fname = fullfile(pth, ['MT0' suffix]);
        create(out.MT0);
    end
    out.MT0.dat(:,:,:) = exp(out.MTw.dat());
end