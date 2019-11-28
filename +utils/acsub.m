function X = acsub(X, mask)
% Extract the autocalibration region of kspace
%
% FORMAT acdata = utils.acsub(kdata)

dim = [size(X) 1];
dim_other = dim(4:end);
dim = dim(1:3);

X = reshape(X, dim(1), dim(2), dim(3), []);

if nargin < 2 || isempty(mask)
    mask  = sum(sum(X, 3),4) > single(eps);
end

hdim = floor(dim/2);
mask_tl = mask(1:hdim(1),1:hdim(2));
mask_br = mask(hdim(1)+1:end,hdim(2)+1:end);
mask_tl = mask_tl(end:-1:1,end:-1:1);

maxm = size(mask_tl,1);
maxn = size(mask_tl,2);
bestm = 0;
bestn = 0;
best  = 0;
m = 1;
n = 1;
while m <= maxm
    while n <= maxn
        if mask_tl(m,n)
            if n*m > best
                bestm = m;
                bestn = n;
                best = n*m;
            end
        else
            maxn = n-1;
        end
        n = n + 1;
    end
    n = 1;
    m = m + 1;
end
minx = hdim(1) - bestm + 1;
miny = hdim(2) - bestn + 1;


maxm = size(mask_br,1);
maxn = size(mask_br,2);
bestm = 0;
bestn = 0;
best  = 0;
m = 1;
n = 1;
while m <= maxm
    while n <= maxn
        if mask_br(m,n)
            if n*m > best
                bestm = m;
                bestn = n;
                best = n*m;
            end
        else
            maxn = n-1;
        end
        n = n + 1;
    end
    n = 1;
    m = m + 1;
end
maxx = hdim(1) + bestm;
maxy = hdim(2) + bestn;

% Extract
idx = {':', ':', ':', ':'};
idx{1} = minx:maxx;
idx{2} = miny:maxy;
S = struct('type', '()', 'subs', {idx});
X = subsref(X, S);

% Reshape
X = reshape(X, [size(X,1) size(X,2) size(X,3) dim_other]);