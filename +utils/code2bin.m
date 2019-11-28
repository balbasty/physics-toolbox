function bin = code2bin(code, length)
% Convert a compact binary mask (an unsigned integer) into a logical mask 
% (a logical vector).
%
% FORMAT bin = utils.code2bin(code, length)
bin = dec2bin(code,length) == '1';
bin = bin(end:-1:1);