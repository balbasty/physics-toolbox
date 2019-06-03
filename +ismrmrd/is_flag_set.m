function is_set = is_flag_set(flags, value, cond)
% FORMAT is_set = is_flag_set(flags, value)
% flags - [uint64] Flag value(s) stored in the hard header
% value - [uint64] Standard value(s) that we hope is in the flag
% cond  - Condition to meet if multiple values are provided: ''/'any'/'all'
%         If empty, return individual true/false for each value
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

if nargin < 3
    cond = '';
end

is_set = de2bi(flags,max(max(value),numel(de2bi(max(flags)))));   % Extract the value of each bit
is_set = is_set(:,value);           % Extract only flags of interest
switch lower(cond)
    case 'any'
        is_set = any(is_set,2);
    case 'all'
        is_set = all(is_set,2);
end

% % This version is probably faster but less easy to read
% is_set = strcmpi(cond, 'all');
% for i=1:numel(value)
%     bitmask = bitshift(uint64(1), uint64(value(i)) - uint64(1));
%     if strcmpi(cond,'all')
%         is_set  = is_set & bsxfun(@bitand, uint64(flags), bitmask) > 0;
%     else
%         is_set  = is_set | bsxfun(@bitand, uint64(flags), bitmask) > 0;
%     end
% end