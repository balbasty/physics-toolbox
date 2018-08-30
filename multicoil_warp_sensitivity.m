function [ws] = multicoil_warp_sensitivity(s, acf)
% FORMAT idx = multicoil_warp_sensitivity(lat, acf)
% lat - [n1 n2 n3] Complete (= unwarped) lattice size
% acf - [f1 f2 f2] Acceleration factor along each dimension
% idx - [n1/f1 n2/f2 n3/f3 n1 n2 n3] Indices to warp the sensitivity
%
% FORMAT ws  = multicoil_warp_sensitivity(s, acf)
% s   - Sensitivity field with lattice [n1 n2 n3 nc]
% acf - [f1 f2 f2] Acceleration factor along each dimension
% ws  - Warped sensitivity field with lattice [n1/f1 n2/f2 n3/f3 n1 n2 n3 nc]

if numel(size(s)) == 2 && size(s,1) == 1
    % Lattice case
    lat = s;
    
else
    % Sensitivity case
end