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
%                        /!\ caipi is experimental
%
% OPTIONS (KEYWORDS)
% ------------------
% layout - Output memory layout: ['box']/'expand'/'compact'
%          * box:     crop input indices with no data
%          * expand:  keep input indices (might use up a lot of memory)
%          * compact: Do not store zeros (indices are not "right")
% order  - Output dimensions ordering
%          [{'ch' 'rd' 'k1' 'k2' 'av' 'sl' 'ct' 'ph' 'rp' 'st' 'sg'}]
%
% OUTPUT
% ------
% dat - A k-space volume.

% -------------------------------------------------------------------------
% Parse input arguments
default_order = {'ch' 'rd' 'k1' 'k2' 'av' 'sl' 'ct' 'ph' 'rp' 'st' 'sg'};
p = inputParser;
p.FunctionName = 'ismrmrd_read';
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
p.addParameter('layout',              'box', @ischar);
p.addParameter('order',       default_order, @iscell);
p.addParameter('subpart',                '', @ischar);
p.parse(varargin{:});
k1_out = p.Results.kspace_encode_step_1;
k2_out = p.Results.kspace_encode_step_2;
av_out = p.Results.average;
sl_out = p.Results.slice;
ct_out = p.Results.contrast;
ph_out = p.Results.phase;
rp_out = p.Results.repetition;
st_out = p.Results.set;
sg_out = p.Results.segment;
rd_out = p.Results.readout;
ch_out = p.Results.channel;
layout  = p.Results.layout;
subpart = p.Results.subpart;
order   = p.Results.order;

flags_remove = {}; % Unused for now, might be useful?
flags_out    = {}; % Unused for now, might be useful?

% -------------------------------------------------------------------------
% Read soft and hard headers
fprintf('Read header\n');
softHeader = ismrmrd_xml(fname);
softHeader = softHeader.ismrmrdHeader;
[hardHeader, nblines] = ismrmrd_read_member(fname, 'head');
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
% Specific subparts
switch lower(subpart)
    case {'autocalib' 'ac'}
        fprintf('Find autocalibration lines\n');
        try
            idx_k1    = hardHeader.idx.kspace_encode_step_1;
            centre_k1 = limits.kspace_encoding_step_1.centre;
            idx_k2    = hardHeader.idx.kspace_encode_step_2;
            centre_k2 = limits.kspace_encoding_step_2.centre;
            
            ack1name       = 'EmbeddedRefLinesE1';
            ack2name       = 'EmbeddedRefLinesE2';
            userparams     = softHeader.userParameters.userParameterLong;
            usernames      = cellfun(@(X) X.name,  userparams, 'UniformOutput', false);
            uservalues     = cellfun(@(X) X.value, userparams);
            ac_nblines_k1  = uservalues(strcmpi(usernames,ack1name));
            ac_nblines_k2  = uservalues(strcmpi(usernames,ack2name));
            
            mask_k1 = idx_k1 >= centre_k1 - ac_nblines_k1/2 ...
                    & idx_k1 <  centre_k1 + ac_nblines_k1/2;
            mask_k2 = idx_k2 >= centre_k2 - ac_nblines_k2/2 ...
                    & idx_k2 <  centre_k2 + ac_nblines_k2/2;
            mask    = mask & mask_k1(:) & mask_k2(:);
            clear mask_k1 mask_k2
        catch
            error('Could not extract autocalibration lines');
        end
        
    case 'cartesian'
        fprintf('Find cartesian lines\n');
        try
            idx_k1  = hardHeader.idx.kspace_encode_step_1;
            min_k1  = limits.kspace_encoding_step_1.minimum;
            max_k1  = limits.kspace_encoding_step_1.maximum;
            idx_k2  = hardHeader.idx.kspace_encode_step_2;
            min_k2  = limits.kspace_encoding_step_2.minimum;
            max_k2  = limits.kspace_encoding_step_2.maximum;
            af_k1   = softHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
            af_k2   = softHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2;
            
            mask_k1 = ismember(idx_k1, min_k1:af_k1:max_k1);
            mask_k2 = ismember(idx_k2, min_k2:af_k2:max_k2);
            mask    = mask & mask_k1(:) & mask_k2(:);
            clear mask_k1 mask_k2 idx_k1 idx_k2 
        catch
            error('Could not extract cartesian lines');
        end
        
    case 'caipi'
        fprintf('Find CAIPI lines\n');
        warning('CAIPI extraction was never tested.')
        try
            idx_k1 = hardHeader.idx.kspace_encode_step_1;
            min_k1 = limits.kspace_encoding_step_1.minimum;
            max_k1 = limits.kspace_encoding_step_1.maximum;
            idx_k2 = hardHeader.idx.kspace_encode_step_2;
            min_k2 = limits.kspace_encoding_step_2.minimum;
            max_k2 = limits.kspace_encoding_step_2.maximum;
            af_k1   = softHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
            af_k2   = softHeader.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2;

            mask_even_k1 = ismember(idx_k1, min_k1:af_k1:max_k1);
            mask_even_k2 = ismember(idx_k2, min_k2:2*af_k2:max_k2);
            mask_odd_k1  = ismember(idx_k1, min_k1+1:af_k1:max_k1);
            mask_odd_k2  = ismember(idx_k2, min_k2+2:2*af_k2:max_k2);
            mask         = mask & ( (mask_even_k1(:) & mask_even_k2(:)) | ...
                                    (mask_odd_k1(:)  & mask_odd_k2(:))  );
            clear mask_even_k1 mask_even_k2 mask_odd_k1 mask_odd_k2 idx_k1 idx_k2 
        catch
            error('Could not extract CAIPI lines');
        end
end

% -------------------------------------------------------------------------
% Select indices
fprintf('Select lines\n');
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
flags = struct( ...
    'ACQ_FIRST_IN_ENCODE_STEP1',                1, ...
    'ACQ_LAST_IN_ENCODE_STEP1',                 2, ...
    'ACQ_FIRST_IN_ENCODE_STEP2',                3, ...
    'ACQ_LAST_IN_ENCODE_STEP2',                 4, ...
    'ACQ_FIRST_IN_AVERAGE',                     5, ...
    'ACQ_LAST_IN_AVERAGE',                      6, ...
    'ACQ_FIRST_IN_SLICE',                       7, ...
    'ACQ_LAST_IN_SLICE',                        8, ...
    'ACQ_FIRST_IN_CONTRAST',                    9, ...
    'ACQ_LAST_IN_CONTRAST',                    10, ...
    'ACQ_FIRST_IN_PHASE',                      11, ...
    'ACQ_LAST_IN_PHASE',                       12, ...
    'ACQ_FIRST_IN_REPETITION',                 13, ...
    'ACQ_LAST_IN_REPETITION',                  14, ...
    'ACQ_FIRST_IN_SET',                        15, ...
    'ACQ_LAST_IN_SET',                         16, ...
    'ACQ_FIRST_IN_SEGMENT',                    17, ...
    'ACQ_LAST_IN_SEGMENT',                     18, ...
    'ACQ_IS_NOISE_MEASUREMENT',                19, ...
    'ACQ_IS_PARALLEL_CALIBRATION',             20, ...
    'ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', 21, ...
    'ACQ_IS_REVERSE',                          22, ...
    'ACQ_IS_NAVIGATION_DATA',                  23, ...
    'ACQ_IS_PHASECORR_DATA',                   24, ...
    'ACQ_LAST_IN_MEASUREMENT',                 25, ...
    'ACQ_IS_HPFEEDBACK_DATA',                  26, ...
    'ACQ_IS_DUMMYSCAN_DATA',                   27, ...
    'ACQ_IS_RTFEEDBACK_DATA',                  28, ...
    'ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA',   29, ...
    'ACQ_USER1',                               57, ...
    'ACQ_USER2',                               58, ...
    'ACQ_USER3',                               59, ...
    'ACQ_USER4',                               60, ...
    'ACQ_USER5',                               61, ...
    'ACQ_USER6',                               62, ...
    'ACQ_USER7',                               63, ...
    'ACQ_USER8',                               64);

function X = flag_to_id(X)
    if ischar(X), X = flags.(X); end
end
if ~isempty(flags_remove)
    if ~iscell(flags_remove), flags_remove = {flags_remove}; end
    flags_remove = cellfun(@flag_to_id, flags_remove);
    mask = mask & reshape(~ismember(hardHeader.flags, flags_remove), [], 1);
end
if ~isempty(flags_out)
    if ~iscell(flags_out), flags_out = {flags_out}; end
    flags_out = cellfun(@flag_to_id, flags_out);
    mask = mask & reshape(ismember(hardHeader.flags, flags_out), [], 1);
end

% -------------------------------------------------------------------------
% Read selected lines from dataset
fprintf('Read selected lines\n');
dataLines  = ismrmrd_read_member(fname, {'data' 'head'}, find(mask));
rd_size    = max(dataLines.head.number_of_samples);
ch_size    = max(dataLines.head.available_channels);
dataLines  = cat(2, dataLines.data{:});
dataLines  = reshape(dataLines, 2, rd_size, ch_size, []);

% -------------------------------------------------------------------------
% Remove unwanted channels and samples 
% + convert to complex
% + orientation [ch rd other]
fprintf('Select samples and channels\n');
if isempty(rd_out), rd_out = 1:rd_size; else, rd_size = numel(rd_out); end
if isempty(ch_out), ch_out = 1:ch_size; else, ch_size = numel(ch_out); end
dataLines = dataLines(:,rd_out,ch_out,:);
dataLines = permute(dataLines(1,:,:,:) + 1i * dataLines(2,:,:,:), [3 2 4 1]);

% -------------------------------------------------------------------------
% Allocate output + populate 
fprintf('Format output\n');
idx_out = cell(1,9); % < [k1 k2 av sl ct ph rp st sg]
switch lower(layout)
    case 'expand'
        idx_out{1} = 1 + hardHeader.idx.kspace_encode_step_1(mask);
        idx_out{2} = 1 + hardHeader.idx.kspace_encode_step_2(mask);
        idx_out{3} = 1 + hardHeader.idx.average(mask);
        idx_out{4} = 1 + hardHeader.idx.slice(mask);
        idx_out{5} = 1 + hardHeader.idx.contrast(mask);
        idx_out{6} = 1 + hardHeader.idx.phase(mask);
        idx_out{7} = 1 + hardHeader.idx.repetition(mask);
        idx_out{8} = 1 + hardHeader.idx.set(mask);
        idx_out{9} = 1 + hardHeader.idx.segment(mask);
    case 'box'
        idx_out{1} = 1 + hardHeader.idx.kspace_encode_step_1(mask);
        idx_out{2} = 1 + hardHeader.idx.kspace_encode_step_2(mask);
        idx_out{3} = 1 + hardHeader.idx.average(mask);
        idx_out{4} = 1 + hardHeader.idx.slice(mask);
        idx_out{5} = 1 + hardHeader.idx.contrast(mask);
        idx_out{6} = 1 + hardHeader.idx.phase(mask);
        idx_out{7} = 1 + hardHeader.idx.repetition(mask);
        idx_out{8} = 1 + hardHeader.idx.set(mask);
        idx_out{9} = 1 + hardHeader.idx.segment(mask);
        
        idx_out{1} = 1 + idx_out{1} - min(idx_out{1});
        idx_out{2} = 1 + idx_out{2} - min(idx_out{2});
        idx_out{3} = 1 + idx_out{3} - min(idx_out{3});
        idx_out{4} = 1 + idx_out{4} - min(idx_out{4});
        idx_out{5} = 1 + idx_out{5} - min(idx_out{5});
        idx_out{6} = 1 + idx_out{6} - min(idx_out{6});
        idx_out{7} = 1 + idx_out{7} - min(idx_out{7});
        idx_out{8} = 1 + idx_out{8} - min(idx_out{8});
        idx_out{9} = 1 + idx_out{9} - min(idx_out{9});
    case 'compact'
        [idx_out{3},~] = find(bsxfun(@eq, hardHeader.idx.average(mask),    unique(hardHeader.idx.average(mask))'));
        [idx_out{4},~] = find(bsxfun(@eq, hardHeader.idx.slice(mask),      unique(hardHeader.idx.slice(mask))'));
        [idx_out{5},~] = find(bsxfun(@eq, hardHeader.idx.contrast(mask),   unique(hardHeader.idx.contrast(mask))'));
        [idx_out{6},~] = find(bsxfun(@eq, hardHeader.idx.phase(mask),      unique(hardHeader.idx.phase(mask))'));
        [idx_out{7},~] = find(bsxfun(@eq, hardHeader.idx.repetition(mask), unique(hardHeader.idx.repetition(mask))'));
        [idx_out{8},~] = find(bsxfun(@eq, hardHeader.idx.set(mask),        unique(hardHeader.idx.set(mask))'));
        [idx_out{9},~] = find(bsxfun(@eq, hardHeader.idx.segment(mask),    unique(hardHeader.idx.segment(mask))'));
        switch lower(subpart)
            case 'caipi'
                error('Compact CAIPI not implemented yet')
            otherwise
                [idx_out{1},~] = find(bsxfun(@eq, hardHeader.idx.kspace_encode_step_1(mask), unique(hardHeader.idx.kspace_encode_step_1(mask))'));
                [idx_out{2},~] = find(bsxfun(@eq, hardHeader.idx.kspace_encode_step_2(mask), unique(hardHeader.idx.kspace_encode_step_2(mask))'));
        end
    otherwise
        error('Unknown layout %s', layout)
end
idx_out = cellfun(@(X) X(:), idx_out, 'UniformOutput', false);
dim_out = [ch_size rd_size cellfun(@max, idx_out)];
dat     = zeros(dim_out, 'like', dataLines);

dat(:,:,sub2ind(dim_out(3:end), idx_out{:})) = dataLines;

% Reorder dimensions
permutation = zeros(1,numel(order));
all_dim     = 1:11;
for i=1:numel(order)
    permutation(i) = find(strcmpi(order{i}, default_order), 1);
end
permutation = [permutation all_dim(~ismember(all_dim, permutation))];
dat         = permute(dat, permutation);

end