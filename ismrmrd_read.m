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
% Parse input arguments
p = inputParser;
p.FunctionName = 'ismrmrd_read';
p.addParameter('channel',               Inf, @isnumeric);
p.addParameter('readout',               Inf, @isnumeric);
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
p.addParameter('maxmem',               5E8, @islogical);
p.parse(varargin{:});
k1 = p.Results.kspace_encode_step_1;
k2 = p.Results.kspace_encode_step_2;
av = p.Results.average;
sl = p.Results.slice;
ct = p.Results.contrast;
ph = p.Results.phase;
rp = p.Results.repetition;
st = p.Results.set;
sg = p.Results.segment;
rd = p.Results.readout;
ch = p.Results.channel;
compact = p.Results.compact;
maxmem  = p.Results.maxmem;

% -------------------------------------------------------------------------
% Read soft header
% fprintf('Read soft header\n');
% softhead = ismrmrd_xml(fname);
info     = h5info(fname, '/dataset/data');
nlines   = info.Dataspace.Size;

% -------------------------------------------------------------------------
% Reading strategy
dat1 = h5read(fname, '/dataset/data', 1, 1);
linelength = numel(dat1.data{1});
clear dat1
linestep = floor(maxmem/(linelength*4));

% -------------------------------------------------------------------------
% Read available indices
fprintf('Read hard header:   0%%');
k1steps     = [];
k2steps     = [];
averages    = [];
slices      = [];
contrasts   = [];
phases      = [];
repetitions = [];
sets        = [];
segments    = [];
pct_prev    = 0;
for i=1:linestep:nlines
    pct = floor(100*i/nlines);
    if pct ~= pct_prev
        pct_prev = pct;
        fprintf('\b\b\b\b%3d%%', pct);
    end
    nread       = min(linestep, nlines-i);
    dat1        = h5read(fname, '/dataset/data', i, nread);
    k1steps     = unique([k1steps       dat1.head.idx.kspace_encode_step_1(:)']);
    k2steps     = unique([k2steps       dat1.head.idx.kspace_encode_step_2(:)']);
    averages    = unique([averages      dat1.head.idx.average(:)']);
    slices      = unique([slices        dat1.head.idx.slice(:)']);
    contrasts   = unique([contrasts     dat1.head.idx.contrast(:)']);
    phases      = unique([phases        dat1.head.idx.phase(:)']);
    repetitions = unique([repetitions   dat1.head.idx.repetition(:)']);
    sets        = unique([sets          dat1.head.idx.set(:)']);
    segments    = unique([segments      dat1.head.idx.segment(:)']);
end
fprintf('\n');
% Convert to matlab 1-indices
k1steps     = k1steps     + 1;
k2steps     = k2steps     + 1;
averages    = averages    + 1;
slices      = slices      + 1;
contrasts   = contrasts   + 1;
phases      = phases      + 1;
repetitions = repetitions + 1;
sets        = sets        + 1;
segments    = segments    + 1;

% -------------------------------------------------------------------------
% Read index range inside each line
dat1      = h5read(fname, '/dataset/data', 1, 1);
nreadouts = dat1.head.number_of_samples;
nchannels = dat1.head.available_channels;
channels  = 1:nchannels;
readouts  = 1:nreadouts;

% -------------------------------------------------------------------------
% Default indices
if ch == Inf, ch = channels;    end
if rd == Inf, rd = readouts;    end
if k1 == Inf, k1 = k1steps;     end
if k2 == Inf, k2 = k2steps;     end
if av == Inf, av = averages;    end
if sl == Inf, sl = slices;      end
if ct == Inf, ct = contrasts;   end
if ph == Inf, ph = phases;      end
if rp == Inf, rp = repetitions; end
if st == Inf, st = sets;        end
if sg == Inf, sg = segments;    end

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
dim = [numel(ch) numel(rd) getdim(k1) getdim(k2) numel(av) ...
       numel(sl) numel(ct) numel(ph)  numel(rp)  numel(st) numel(sg)];
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
    
    includedind = {k1, k2, av, sl, ct, ph, rp, st, sg};
    foundind    = {k11, k21, av1, sl1, ct1, ph1, rp1, st1, sg1};
    isincluded  = @(x,y) any(bsxfun(@eq, reshape(x, [], 1), ...
                                         reshape(y, 1, [])), 2);
    msk = cellfun(isincluded, foundind, includedind, 'UniformOutput', false);
    msk = all([msk{:}],2);
    
    line = cat(2, dat1.data{msk});
    line = reshape(line, [2 nreadouts nchannels sum(msk)]);
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
        
     dat(:,:,linind) = line(ch,rd,:);
                   
end
fprintf('\n');

end