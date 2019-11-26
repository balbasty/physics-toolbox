function bin = code2bin(code, length)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.code2bin
%--------------------------------------------------------------------------
% FORMAT bin = gmm.lib.code2bin(code, length)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

bin = dec2bin(code,length) == '1';
bin = bin(end:-1:1);