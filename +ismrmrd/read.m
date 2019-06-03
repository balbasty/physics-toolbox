function dat = read(fname, varargin)
% Read (partially) an ismrmrd dataset
%
% FORMAT dat = ismrmrd.read(fname, ...)
%
% INDEX SELECTION (KEYWORDS)
% --------------------------
% channel              - Indices of channels to read
% readout              - Indices of readout samples to read
% kspace_encode_step_1 - Indices of first phase encoding dir to read
% kspace_encode_step_2 - Indices of second phase encoding dir to read
% average              - Indices of averages to read
% slice                - Indices of slices to read
% contrast             - Indices of contrasts/echoes to read
% phase                - Indices of phases to read
% repetition           - Indices of repetitions to read
% set                  - Indices of sets to read
% segment              - Indices of segments to read
% 
% FLAG SELECTION (KEYWORDS)
% -------------------------
% subpart   - Specific cases
%             ['kspace']/'imaging'/'calibration'/'noise'/'navigation'/'perso'
% flags_in  - Selected lines must have one of these flags:
%             * 'kspace':     {} [unused]
%             * 'imaging':    {} [unused]
%             * 'calib':      {20 21}
%             * 'noise':      {19}
%             * 'navigation': {23 24}
%             * 'perso':      {} [unused]
% flags_out - Selected lines cannot have one of these flags 
%             * 'kspace':     {19 23 24 26 27 28 29}
%             * 'imaging':    {19 20 23 24 26 27 28 29}
%             * 'calib':      {19 23 24 26 27 28 29}
%             * 'noise':      {}
%             * 'navigation': {}
%             * 'perso':      {}
%             
% See ismrmrd.acquisition_flags to get the list of available flags.
% flags_in/flags_out must be cells of strings or integers.
%
% OPTIONS (KEYWORDS)
% ------------------
% layout  - Output memory layout: 'acq'/'box'/'tight'/'vector'
%           * 'kspace':     'acq'       (use acquisition matrix)
%           * 'imaging':    'acq'       (use acquisition matrix)
%           * 'calib':      'box'       (crop at largest absolute freq)
%           * 'noise':      'vector'    (do not organise the data)
%           * 'navigation': 'vector'    (do not organise the data)
%           * 'perso':      'vector'    (do not organise the data)
% order   - Output dimensions ordering
%           [{'ch' 'rd' 'k1' 'k2' 'av' 'sl' 'ct' 'ph' 'rp' 'st' 'sg'}]
% verbose - Speak while reading [false]
%
% OUTPUT
% ------
% dat  - A k-space volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse input arguments
default_order = {'ch' 'rd' 'k1' 'k2' 'av' 'sl' 'ct' 'ph' 'rp' 'st' 'sg'};
p = inputParser;
p.FunctionName = 'ismrmrd.read';
p.addParameter('channel',               [],  @isnumeric);
p.addParameter('readout',               [],  @(X) isnumeric(X) || ischar(X));
p.addParameter('kspace_encode_step_1',  [],  @isnumeric);
p.addParameter('kspace_encode_step_2',  [],  @isnumeric);
p.addParameter('average',               [],  @isnumeric);
p.addParameter('slice',                 [],  @isnumeric);
p.addParameter('contrast',              [],  @isnumeric);
p.addParameter('phase',                 [],  @isnumeric);
p.addParameter('repetition',            [],  @isnumeric);
p.addParameter('set',                   [],  @isnumeric);
p.addParameter('segment',               [],  @isnumeric);
p.addParameter('subpart',          'kspace', @ischar);
p.addParameter('flags_in',               {}, @(X) iscell(X) || ischar(X) || isscalar(X));
p.addParameter('flags_out',              {}, @(X) iscell(X) || ischar(X) || isscalar(X));
p.addParameter('layout',                 '', @ischar);
p.addParameter('order',       default_order, @iscell);
p.addParameter('verbose',                 0, @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.parse(varargin{:});
k1_out    = p.Results.kspace_encode_step_1;
k2_out    = p.Results.kspace_encode_step_2;
av_out    = p.Results.average;
sl_out    = p.Results.slice;
ct_out    = p.Results.contrast;
ph_out    = p.Results.phase;
rp_out    = p.Results.repetition;
st_out    = p.Results.set;
sg_out    = p.Results.segment;
rd_out    = p.Results.readout;
ch_out    = p.Results.channel;
layout    = p.Results.layout;
subpart   = p.Results.subpart;
flags_in  = p.Results.flags_in;
flags_out = p.Results.flags_out;
order     = p.Results.order;
verbose   = p.Results.verbose;


% -------------------------------------------------------------------------
% Preprocess options
if ~iscell(flags_in)
    flags_in = {flags_in};
end
switch lower(subpart)
    case {'kspace' 'imaging' 'perso'}
    case 'calibration'
        flags_in = [flags_in {20 21}];
    case 'noise'
        flags_in = [flags_in {19}];
    case 'navigation'
        flags_in = [flags_in {23 24}];
    otherwise
        warning('Unknown subpart %s', subpart);
end
if ~iscell(flags_out)
    flags_out = {flags_out};
end
switch lower(subpart)
    case {'kspace' 'calibration'}
        flags_out = [flags_out {19 23 24 26 27 28 29}];
    case 'imaging'
        flags_out = [flags_out {19 20 23 24 26 27 28 29}];
    case {'noise' 'navigation' 'perso'}
    otherwise
        warning('Unknown subpart %s', subpart);
end
if isempty(layout)
    switch lower(subpart)
        case {'kspace' 'imaging'}
            layout = 'acq';
        case 'calibration'
            layout = 'box';
        case {'noise' 'navigation' 'perso'}
            layout = 'vector';
        otherwise
            warning('Unknown subpart %s', subpart);
    end
end

% -------------------------------------------------------------------------
% Read soft and hard headers
if verbose, fprintf('Read header\n'); end
softHeader = ismrmrd.xml(fname);
softHeader = softHeader.ismrmrdHeader;
[hardHeader, nblines] = ismrmrd.read_member(fname, 'head');
hardHeader = hardHeader.head;

% -------------------------------------------------------------------------
% Store limits along dimensions
limits                = softHeader.encoding.encodingLimits;
limits.samples.min    = 0;
limits.samples.max    = max(hardHeader.number_of_samples) - 1;
limits.samples.center = limits.samples.max/2;
limits.channels.min   = 1;
limits.channels.max   = max(hardHeader.available_channels);

% -------------------------------------------------------------------------
% Prepare mask of lines to read
mask = ones(nblines, 1, 'logical');

% -------------------------------------------------------------------------
% Select indices
if verbose, fprintf('Select lines\n'); end
if ~isempty(k1_out)
    mask = mask & reshape(ismember(hardHeader.idx.kspace_encode_step_1, k1_out-1), [], 1);
end
if ~isempty(k2_out)
    mask = mask & reshape(ismember(hardHeader.idx.kspace_encode_step_2, k2_out-1), [], 1);
end
if ~isempty(av_out)
    mask = mask & reshape(ismember(hardHeader.idx.average, av_out-1), [], 1);
end
if ~isempty(sl_out)
    mask = mask & reshape(ismember(hardHeader.idx.slice, sl_out-1), [], 1);
end
if ~isempty(ct_out)
    mask = mask & reshape(ismember(hardHeader.idx.contrast, ct_out-1), [], 1);
end
if ~isempty(ph_out)
    mask = mask & reshape(ismember(hardHeader.idx.phase, ph_out-1), [], 1);
end
if ~isempty(rp_out)
    mask = mask & reshape(ismember(hardHeader.idx.repetition, rp_out-1), [], 1);
end
if ~isempty(st_out)
    mask = mask & reshape(ismember(hardHeader.idx.set, st_out-1), [], 1);
end
if ~isempty(sg_out)
    mask = mask & reshape(ismember(hardHeader.idx.segment, sg_out-1), [], 1);
end

% -------------------------------------------------------------------------
% Select flags
flags = ismrmrd.acquisition_flags;

function X = flag_to_id(X)
    if ischar(X), X = flags.(X); end
end
if ~isempty(flags_out)
    if ~iscell(flags_out), flags_out = {flags_out}; end
    flags_out = cellfun(@flag_to_id, flags_out);
    mask = mask & reshape(~ismrmrd.is_flag_set(hardHeader.flags, flags_out, 'any'), [], 1);
end
if ~isempty(flags_in)
    if ~iscell(flags_in), flags_in = {flags_in}; end
    flags_in = cellfun(@flag_to_id, flags_in);
    mask = mask & reshape(ismrmrd.is_flag_set(hardHeader.flags, flags_in, 'any'), [], 1);
end

% -------------------------------------------------------------------------
% Read selected lines from dataset
if verbose, fprintf('Read selected lines\n'); end
dataLines  = ismrmrd.read_member(fname, {'data' 'head'}, find(mask)-1);
rd_size    = max(dataLines.head.number_of_samples);
ch_size    = max(dataLines.head.available_channels);
dataLines  = cat(2, dataLines.data{:});
dataLines  = reshape(dataLines, 2, rd_size, ch_size, []); % [r/i rd ch ...]
dataLines  = permute(dataLines, [3 2 4 1]);               % [ch rd ... r/i]

% -------------------------------------------------------------------------
% Remove unwanted channels and samples 
% + convert to complex
% + orientation [ch rd other]
if verbose, fprintf('Select samples and channels\n'); end
if isempty(rd_out), rd_out = 1:rd_size; else, rd_size = numel(rd_out); end
if isempty(ch_out), ch_out = 1:ch_size; else, ch_size = numel(ch_out); end
if numel(ch_out) ~= size(dataLines,1) || numel(rd_out) ~= size(dataLines,2)
    dataLines = dataLines(ch_out,rd_out,:,:);
end
dataLines = complex(dataLines(:,:,:,1),dataLines(:,:,:,2));


% -------------------------------------------------------------------------
% Allocate output + populate 
if verbose, fprintf('Format output\n'); end
idx_out = cell(1,9); % < [k1 k2 av sl ct ph rp st sg]

% Convert 0-indexing into 1-indexing
idx_out{1} = 1 + int32(hardHeader.idx.kspace_encode_step_1(mask));
idx_out{2} = 1 + int32(hardHeader.idx.kspace_encode_step_2(mask));
idx_out{3} = 1 + int32(hardHeader.idx.average(mask));
idx_out{4} = 1 + int32(hardHeader.idx.slice(mask));
idx_out{5} = 1 + int32(hardHeader.idx.contrast(mask));
idx_out{6} = 1 + int32(hardHeader.idx.phase(mask));
idx_out{7} = 1 + int32(hardHeader.idx.repetition(mask));
idx_out{8} = 1 + int32(hardHeader.idx.set(mask));
idx_out{9} = 1 + int32(hardHeader.idx.segment(mask));

% average/slice/contrast/phase/repetition/set/segment
% -> for these, we do not want to allocate indices that are not read.
if ~isempty(av_out)
    [idx_out{3},~] = find(bsxfun(@eq, idx_out{3}, av_out(:)')');
    idx_out{3} = reshape(int32(idx_out{3}), [], 1);
end
if ~isempty(sl_out)
    [idx_out{4},~] = find(bsxfun(@eq, idx_out{4}, sl_out(:)')');
    idx_out{4} = reshape(int32(idx_out{4}), [], 1);
end
if ~isempty(ct_out)
    [idx_out{5},~] = find(bsxfun(@eq, idx_out{5}, ct_out(:)')');
    idx_out{5} = reshape(int32(idx_out{5}), [], 1);
end
if ~isempty(ph_out)
    [idx_out{6},~] = find(bsxfun(@eq, idx_out{6}, ph_out(:)')');
    idx_out{6} = reshape(int32(idx_out{6}), [], 1);
end
if ~isempty(rp_out)
    [idx_out{7},~] = find(bsxfun(@eq, idx_out{7}, rp_out(:)')');
    idx_out{7} = reshape(int32(idx_out{7}), [], 1);
end
if ~isempty(st_out)
    [idx_out{8},~] = find(bsxfun(@eq, idx_out{8}, st_out(:)')');
    idx_out{8} = reshape(int32(idx_out{8}), [], 1);
end
if ~isempty(sg_out)
    [idx_out{9},~] = find(bsxfun(@eq, idx_out{9}, sg_out(:)')');
    idx_out{9} = reshape(int32(idx_out{9}), [], 1);
end

% k1/k2
% -> for these, it depends on the mode
switch lower(layout)
    case 'vector'
        % Here, we don't relate lines with k-space locations. 
        % We keep them ordered as they are within each subvolume
        % - k1 is used to number lines within each (numbered) subvolume
        % - k2 is always 1
        idx_out{2}(:) = 1;
        idx_sub = cellfun(@(X) X(:), idx_out(3:end), 'UniformOutput', false);
        dim_sub = cellfun(@max, idx_sub);
        for av=1:dim_sub(1)
        for sl=1:dim_sub(2)
        for ct=1:dim_sub(3)
        for ph=1:dim_sub(4)
        for rp=1:dim_sub(5)
        for st=1:dim_sub(6)
        for sg=1:dim_sub(7)
            sub_mask = (idx_sub{1} == av) ...
                    & (idx_sub{2} == sl) ...
                    & (idx_sub{3} == ct) ...
                    & (idx_sub{4} == ph) ...
                    & (idx_sub{5} == rp) ...
                    & (idx_sub{6} == st) ...
                    & (idx_sub{7} == sg);
            idx_out{1}(sub_mask) = 1:numel(sub_mask);
        end; end; end; end; end; end; end
        k1_size = max(idx_out{1});
        k2_size = 1;
    case 'acq'
        % We use stored indices
        k1_size = softHeader.encoding.encodedSpace.matrixSize.y;
        k2_size = softHeader.encoding.encodedSpace.matrixSize.z;
    case 'box'
        % The min(index) or (centre_index-max(index))
        nb_lower_k1 = int32(limits.kspace_encoding_step_1.center) ...
            - min(idx_out{1}) + 1;
        nb_upper_k1 = max(idx_out{1}) - 1 ...
            - limits.kspace_encoding_step_1.center;
        if nb_lower_k1 >= nb_upper_k1
            shift_k1 = 1 - min(idx_out{1});
        else
            shift_k1 = 1 + nb_upper_k1 - nb_lower_k1 - min(idx_out{1});
        end
        idx_out{1} = idx_out{1} + shift_k1;
        k1_size = max(nb_upper_k1,nb_lower_k1)*2;
        if nb_lower_k1 < nb_upper_k1
            k1_size = k1_size + 1;
        end

        nb_lower_k2 = limits.kspace_encoding_step_2.center ...
            - min(idx_out{2}) + 1;
        nb_upper_k2 = max(idx_out{2}) - 1 ...
            - limits.kspace_encoding_step_2.center;
        if nb_lower_k2 >= nb_upper_k2
            shift_k2 = 1 - min(idx_out{2});
        else
            shift_k2 = 1 + nb_upper_k2 - nb_lower_k2 - min(idx_out{2});
        end
        idx_out{2} = idx_out{2} + shift_k2;
        k2_size = max(nb_upper_k2,nb_lower_k2)*2;
        if nb_lower_k2 < nb_lower_k2
            k2_size = k2_size + 1;
        end
    case 'tight'
        % The smallest read index becomes 1
        idx_out{1} = 1 + idx_out{1} - min(idx_out{1});
        idx_out{2} = 1 + idx_out{2} - min(idx_out{2});
        k1_size = max(idx_out{1});
        k2_size = max(idx_out{2});
    otherwise
        error('Unknown layout %s', layout)
end
idx_out = cellfun(@(X) X(:), idx_out, 'UniformOutput', false);
dim_out = [ch_size rd_size k1_size k2_size cellfun(@max, idx_out(3:end))];
if strcmpi(layout, 'acq')
    dim_out(2) = softHeader.encoding.encodedSpace.matrixSize.x;
    dim_out(3) = softHeader.encoding.encodedSpace.matrixSize.y;
    if unique(hardHeader.idx.slice) > 1
        dim_out(6) = softHeader.encoding.encodedSpace.matrixSize.z;
    else
        dim_out(4) = softHeader.encoding.encodedSpace.matrixSize.z;
    end
end
dat     = zeros(dim_out, 'like', dataLines);

dat(:,:,sub2ind(dim_out(3:end), idx_out{:})) = dataLines;


% -------------------------------------------------------------------------
% Reorder dimensions
permutation = zeros(1,numel(order));
all_dim     = 1:11;
for i=1:numel(order)
    permutation(i) = find(strcmpi(order{i}, default_order), 1);
end
permutation = [permutation all_dim(~ismember(all_dim, permutation))];
dat         = permute(dat, permutation);

end