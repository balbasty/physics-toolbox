function s = model(dim, vs, ncoils)

s      = struct;
s.mean = zeros(dim, 'like', single(1i));
s.sens = zeros([dim ncoils], 'like', single(1i));
s.b1   = zeros(dim, 'single');
s.vs   = vs;