function dat = ismrmrd_read(fname, varargin)
% Read (partially) an ismrmrd dataset
%
% FORMAT dat = ismrmrd_read(fname, ...)
%
% INDEX SELECTION (KEYWORDS)
% --------------------------
% channel              - Indices of channels to read
% readout              - Indices of readout samples to read
% kspace_encode_step_1 - Indices of first phase encoding dir 
% kspace_encode_step_2 - Indices of second phase encoding dir 
% average              - Indices of averages to read
% slice                - Indices of slices to read
% contrast             - Indices of contrasts/echoes to read
% phase                - Indices of phases to read
% repetition           - Indices of repetitions to read
% set                  - Indices of sets to read
% segment              - Indices of segments to read
% subpart              - Specific cases ('autocalib'/'cartesian'/'caipi')
%
% OPTIONS (KEYWORDS)
% ------------------
% compact - Return a compact volume (= do not zero-fill)    [false]
% maxmem  - Peak amount of RAM to use when reading          [500M]
%
% OUTPUT
% ------
% dat - A k-space volume. Dimensions (from most rapidly changing) are:
%           1: channels
%           2: readout samples
%           3: first phase encoding direction
%           4: second phase encoding direction
%           5: averages
%           6: slices
%           7: contrasts or echoes
%           8: repetitions
%           9: sets
%          10: segments

% -------------------------------------------------------------------------
% Prepare header information
head = struct;

% -------------------------------------------------------------------------
% Parse input arguments
p = inputParser;
p.FunctionName = 'ismrmrd_read';
p.addParameter('channel',               Inf, @isnumeric);
p.addParameter('readout',               Inf, @(X) isnumeric(X) || ischar(X));
p.addParameter('kspace_encode_step_1',  Inf, @isnumeric);
p.addParameter('kspace_encode_step_2',  Inf, @isnumeric);
p.addParameter('average',               Inf, @isnumeric);
p.addParameter('slice',                 Inf, @isnumeric);
p.addParameter('contrast',              Inf, @isnumeric);
p.addParameter('phase',                 Inf, @isnumeric);
p.addParameter('repetition',            Inf, @isnumeric);
p.addParameter('set',                   Inf, @isnumeric);
p.addParameter('segment',               Inf, @isnumeric);
p.addParameter('compact',             false, @islogical);
p.addParameter('maxmem',                5E8, @islogical);
p.addParameter('subpart',                '', @ischar);
p.parse(varargin{:});
head.k1.out = p.Results.kspace_encode_step_1;
head.k2.out = p.Results.kspace_encode_step_2;
head.av.out = p.Results.average;
head.sl.out = p.Results.slice;
head.ct.out = p.Results.contrast;
head.ph.out = p.Results.phase;
head.rp.out = p.Results.repetition;
head.st.out = p.Results.set;
head.sg.out = p.Results.segment;
head.rd.out = p.Results.readout;
head.ch.out = p.Results.channel;
compact     = p.Results.compact;
maxmem      = p.Results.maxmem;
subpart     = p.Results.subpart;

% -------------------------------------------------------------------------
% Read soft header
infoh5   = h5info(fname, '/dataset/data');
nlines   = infoh5.Dataspace.Size;
fprintf('Read soft header\n');
softhead = ismrmrd_xml(fname);

head_limits = softhead.ismrmrdHeader.encoding.encodingLimits;

head.k1.lim.min    = 1 + str2double(head_limits.kspace_encoding_step_1.minimum.get_data);
head.k1.lim.max    = 1 + str2double(head_limits.kspace_encoding_step_1.maximum.get_data);
head.k1.lim.center = 1 + str2double(head_limits.kspace_encoding_step_1.center.get_data);
head.k2.lim.min    = 1 + str2double(head_limits.kspace_encoding_step_2.minimum.get_data);
head.k2.lim.max    = 1 + str2double(head_limits.kspace_encoding_step_2.maximum.get_data);
head.k2.lim.center = 1 + str2double(head_limits.kspace_encoding_step_2.center.get_data);
head.av.lim.min    = 1 + str2double(head_limits.average.minimum.get_data);
head.av.lim.max    = 1 + str2double(head_limits.average.maximum.get_data);
head.av.lim.center = 1 + str2double(head_limits.average.center.get_data);
head.sl.lim.min    = 1 + str2double(head_limits.slice.minimum.get_data);
head.sl.lim.max    = 1 + str2double(head_limits.slice.maximum.get_data);
head.sl.lim.center = 1 + str2double(head_limits.slice.center.get_data);
head.ct.lim.min    = 1 + str2double(head_limits.contrast.minimum.get_data);
head.ct.lim.max    = 1 + str2double(head_limits.contrast.maximum.get_data);
head.ct.lim.center = 1 + str2double(head_limits.contrast.center.get_data);
head.ph.lim.min    = 1 + str2double(head_limits.phase.minimum.get_data);
head.ph.lim.max    = 1 + str2double(head_limits.phase.maximum.get_data);
head.ph.lim.center = 1 + str2double(head_limits.phase.center.get_data);
head.rp.lim.min    = 1 + str2double(head_limits.repetition.minimum.get_data);
head.rp.lim.max    = 1 + str2double(head_limits.repetition.maximum.get_data);
head.rp.lim.center = 1 + str2double(head_limits.repetition.center.get_data);
head.st.lim.min    = 1 + str2double(head_limits.set.minimum.get_data);
head.st.lim.max    = 1 + str2double(head_limits.set.maximum.get_data);
head.st.lim.center = 1 + str2double(head_limits.set.center.get_data);
head.sg.lim.min    = 1 + str2double(head_limits.segment.minimum.get_data);
head.sg.lim.max    = 1 + str2double(head_limits.segment.maximum.get_data);
head.sg.lim.center = 1 + str2double(head_limits.segment.center.get_data);

% -------------------------------------------------------------------------
% Reading strategy
dat1 = h5read(fname, '/dataset/data', 1, 1);
linelength = numel(dat1.data{1});
clear dat1
linestep = floor(maxmem/(linelength*4));

% -------------------------------------------------------------------------
% Read available indices
fprintf('Read hard header:   0%%');
head.k1.smp = [];
head.k2.smp = [];
head.av.smp = [];
head.sl.smp = [];
head.ct.smp = [];
head.ph.smp = [];
head.rp.smp = [];
head.st.smp = [];
head.sg.smp = [];
pct_prev    = 0;
for i=1:linestep:nlines
    pct = floor(100*i/nlines);
    if pct ~= pct_prev
        pct_prev = pct;
        fprintf('\b\b\b\b%3d%%', pct);
    end
    nread       = min(linestep, nlines-i);
    dat1        = h5read(fname, '/dataset/data', i, nread);
    head.k1.smp = unique([head.k1.smp dat1.head.idx.kspace_encode_step_1(:)']);
    head.k2.smp = unique([head.k2.smp dat1.head.idx.kspace_encode_step_2(:)']);
    head.av.smp = unique([head.av.smp dat1.head.idx.average(:)']);
    head.sl.smp = unique([head.sl.smp dat1.head.idx.slice(:)']);
    head.ct.smp = unique([head.ct.smp dat1.head.idx.contrast(:)']);
    head.ph.smp = unique([head.ph.smp dat1.head.idx.phase(:)']);
    head.rp.smp = unique([head.rp.smp dat1.head.idx.repetition(:)']);
    head.st.smp = unique([head.st.smp dat1.head.idx.set(:)']);
    head.sg.smp = unique([head.sg.smp dat1.head.idx.segment(:)']);
end
fprintf('\n');
% Convert to matlab 1-indices
head.k1.smp = head.k1.smp + 1;
head.k2.smp = head.k2.smp + 1;
head.av.smp = head.av.smp + 1;
head.sl.smp = head.sl.smp + 1;
head.ct.smp = head.ct.smp + 1;
head.ph.smp = head.ph.smp + 1;
head.rp.smp = head.rp.smp + 1;
head.st.smp = head.st.smp + 1;
head.sg.smp = head.sg.smp + 1;

% -------------------------------------------------------------------------
% Read index range inside each line
dat1            = h5read(fname, '/dataset/data', 1, 1);
head.rd.nb      = dat1.head.number_of_samples;
head.rd.lim.min = 1;
head.rd.lim.max = head.rd.nb;
head.rd.smp     = head.rd.lim.min:head.rd.lim.max;
head.ch.nb      = dat1.head.available_channels;
head.ch.lim.min = 1;
head.ch.lim.max = head.ch.nb;
head.ch.smp     = head.ch.lim.min:head.ch.lim.max;

% -------------------------------------------------------------------------
% Specific subparts
switch subpart
    case 'autocalib'
        try
            ack1name       = 'EmbeddedRefLinesE1';
            ack2name       = 'EmbeddedRefLinesE2';
            userparams     = {softhead.ismrmrdHeader.userParameters.userParameterLong.name};
            userparams     = cellfun(@(x) char(x.get_data), userparams, 'UniformOutput', false);
            uservalues     = {softhead.ismrmrdHeader.userParameters.userParameterLong.value};
            head.k1.ac.nb  = str2double(uservalues{strcmpi(userparams,ack1name)}.get_data);
            head.k2.ac.nb  = str2double(uservalues{strcmpi(userparams,ack2name)}.get_data);
            head.k1.ac.min = head.k1.lim.center - head.k1.ac.nb/2;
            head.k1.ac.max = head.k1.lim.center + head.k1.ac.nb/2 - 1;
            head.k2.ac.min = head.k2.lim.center - head.k2.ac.nb/2;
            head.k2.ac.max = head.k2.lim.center + head.k2.ac.nb/2 - 1;
            head.k1.out    = head.k1.ac.min:head.k1.ac.max;
            head.k2.out    = head.k2.ac.min:head.k2.ac.max;
        catch
        end
    case 'cartesian'
        try
            head.k1.af  = str2double(softhead.ismrmrdHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1.get_data);
            head.k2.af  = str2double(softhead.ismrmrdHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2.get_data);
            head.k1.out = head.k1.lim.min:head.k1.af:head.k1.lim.max;
            head.k2.out = head.k2.lim.min:head.k2.af:head.k2.lim.max;
        catch
        end
end

% -------------------------------------------------------------------------
% Default indices
if ischar(head.rd.out)
    if strcmpi(head.rd.out(1), '/')
        sub         = str2double(head.rd.out(2:end));
        len         = head.rd.lim.max - head.rd.lim.min + 1;
        keep        = ceil(len/sub);
        skip        = floor((len-keep)/2);
        head.rd.out = (head.rd.lim.min + skip):(head.rd.lim.min + skip + keep - 1);
    else
        head.rd.out = Inf;
    end
end
if head.ch.out == Inf, head.ch.out = head.ch.smp; end
if head.rd.out == Inf, head.rd.out = head.rd.smp; end
if head.k1.out == Inf, head.k1.out = head.k1.smp; end
if head.k2.out == Inf, head.k2.out = head.k2.smp; end
if head.av.out == Inf, head.av.out = head.av.smp; end
if head.sl.out == Inf, head.sl.out = head.sl.smp; end
if head.ct.out == Inf, head.ct.out = head.ct.smp; end
if head.ph.out == Inf, head.ph.out = head.ph.smp; end
if head.rp.out == Inf, head.rp.out = head.rp.smp; end
if head.st.out == Inf, head.st.out = head.st.smp; end
if head.sg.out == Inf, head.sg.out = head.sg.smp; end

% -------------------------------------------------------------------------
% Allocate output
compactdim = @(x) numel(x);
expanddim  = @(x) x(end)-x(1)+1;
function i = full2comp(x,y)
    [i,~] = find(bsxfun(@eq, x, y)');
    i = reshape(i, [], 1);
end
compactind = @full2comp;
expandind  = @(x,y) reshape(x - y(1) + 1, [], 1);
if compact, getdim = compactdim; getind = compactind;
else,       getdim = expanddim;  getind = expandind;  end
dim = [numel(head.ch.out)  numel(head.rd.out) ...
       getdim(head.k1.out) getdim(head.k2.out) ...
       numel(head.av.out)  numel(head.sl.out) ...
       numel(head.ct.out)  numel(head.ph.out) ...
       numel(head.rp.out)  numel(head.st.out) ...
       numel(head.sg.out)];
dat = zeros(dim, 'like', single(1i));

% -------------------------------------------------------------------------
% Populate output volume
fprintf('Read data:   0%%');
for i=1:linestep:nlines
    pct = floor(100*i/nlines);
    if pct ~= pct_prev
        pct_prev = pct;
        fprintf('\b\b\b\b%3d%%', pct);
    end
    nread       = min(linestep, nlines-i);
    dat1 = h5read(fname, '/dataset/data', i, nread);
    k11  = 1 + dat1.head.idx.kspace_encode_step_1(:);
    k21  = 1 + dat1.head.idx.kspace_encode_step_2(:);
    av1  = 1 + dat1.head.idx.average(:);
    sl1  = 1 + dat1.head.idx.slice(:);
    ct1  = 1 + dat1.head.idx.contrast(:);
    ph1  = 1 + dat1.head.idx.phase(:);
    rp1  = 1 + dat1.head.idx.repetition(:);
    st1  = 1 + dat1.head.idx.set(:);
    sg1  = 1 + dat1.head.idx.segment(:);
    
    includedind = {head.k1.out, head.k2.out, head.av.out, head.sl.out, head.ct.out, head.ph.out, head.rp.out, head.st.out, head.sg.out};
    foundind    = {k11, k21, av1, sl1, ct1, ph1, rp1, st1, sg1};
    isincluded  = @(x,y) any(bsxfun(@eq, reshape(x, [], 1), ...
                                         reshape(y, 1, [])), 2);
    msk = cellfun(isincluded, foundind, includedind, 'UniformOutput', false);
    msk = all([msk{:}],2);
    
    line = cat(2, dat1.data{msk});
    line = reshape(line, [2 head.rd.nb head.ch.nb sum(msk)]);
    line = permute(line(1,:,:,:) + 1i*line(2,:,:,:), [3 2 4 1]);
    
    foundind = [k11, k21, av1, sl1, ct1, ph1, rp1, st1, sg1];
    foundind = foundind(msk,:);
    linind   = sub2ind(dim(3:end), ...
                           getind(foundind(:,1),includedind{1}), ...
                           getind(foundind(:,2),includedind{2}), ...
                       compactind(foundind(:,3),includedind{3}), ...
                       compactind(foundind(:,4),includedind{4}), ...
                       compactind(foundind(:,5),includedind{5}), ...
                       compactind(foundind(:,6),includedind{6}), ...
                       compactind(foundind(:,7),includedind{7}), ...
                       compactind(foundind(:,8),includedind{8}), ...
                       compactind(foundind(:,9),includedind{9}));
        
     dat(:,:,linind) = line(head.ch.out,head.rd.out,:);
                   
end
fprintf('\n');

end