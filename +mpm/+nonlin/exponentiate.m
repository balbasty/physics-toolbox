function out = exponentiate(out,suffix,pth0)
% Exponentiate log-parameters.
%
% FORMAT out = mpm.nonlin.exponentiate(out,[suffix],[path])
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

if isfield(out, 'logA')
    out.A = out.logA;
    if isa(out.A.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.A.dat.fname);
        end
        out.A.dat.fname = fullfile(pth, ['A' suffix]);
        create(out.A);
    end
    out.A.dat(:,:,:) = exp(out.logA.dat());
end
if isfield(out, 'logR1')
    out.R1 = out.logR1;
    if isa(out.R1.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.R1.dat.fname);
        end
        out.R1.dat.fname = fullfile(pth, ['R1' suffix]);
        create(out.R1);
    end
    out.R1.dat(:,:,:) = exp(out.logR1.dat());
end
if isfield(out, 'logR2s')
    out.R2s = out.logR1;
    if isa(out.R2s.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.R2s.dat.fname);
        end
        out.R2s.dat.fname = fullfile(pth, ['R2s' suffix]);
        create(out.R2s);
    end
    out.R2s.dat(:,:,:) = exp(out.logR2s.dat());
end
if isfield(out, 'logMT')
    out.MT = out.logMT;
    if isa(out.MT.dat, 'file_array')
        pth = pth0;
        if isempty(pth0)
            pth = fileparts(out.MT.dat.fname);
        end
        out.MT.dat.fname = fullfile(pth, ['MT' suffix]);
        create(out.MT);
    end
    out.MT.dat(:,:,:) = exp(out.logMT.dat());
end