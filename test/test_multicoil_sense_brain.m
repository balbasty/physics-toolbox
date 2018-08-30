% -------------------------------------------------------------------------
% Parameters
vs    = [8 8 8];    % > Voxel size
itmax = 150;        % > Maximum number of optim iterations
prm   = [0 1E1 0];    % > Mean image regularisation 
                    %   ['abs value' 'membrane E' 'bending E']
tol   = 1E-4;       % > Tolerance: if gain < tol, stop optimistion
verbose = 2;        

% -------------------------------------------------------------------------
% Phantom data
path = fileparts(which('test_multicoil_brain'));
load(fullfile(path,'brain','brain_results_180828.mat'));

% -------------------------------------------------------------------------
% Sampling scheme
dim = size(rho);
af = [4 2 1]; % acceleration factor
[~,dir] = find(af > 1);
msk = zeros(dim(1),dim(2),dim(3),'logical');
msk(1:af(1):end,1:af(2):end,1:af(3):end) = 1;
ctr = [8 8 dim(3)]; % fully sampled centre
if ~isempty(ctr)
    ctrmsk = zeros(size(msk), 'logical');
    ctrmsk(ceil((dim(1)-ctr(1))/2+1):ceil((dim(1)+ctr(1))/2), ...
           ceil((dim(2)-ctr(2))/2+1):ceil((dim(2)+ctr(2))/2), ...
           ceil((dim(3)-ctr(3))/2+1):ceil((dim(3)+ctr(3))/2)) = 1;
    msk(ctrmsk) = 1;
end

% -------------------------------------------------------------------------
% Wrap observed points
x = multicoil_pullwrap(x,msk,dir);

% -------------------------------------------------------------------------
% Starting estimate
rho  = rho(:,:,:,:,1) + 1i * rho(:,:,:,:,2);
rho0 = rho;
if ~isempty(ctr) % If fully sampled centre, start from there
    for d=dir(:)'
        rho = fft(rho, [], d);
        rho = fftshift(rho, d);
    end
    rho(~ctrmsk) = 0;
    for d=dir(:)'
        rho = fftshift(rho, d);
        rho = ifft(rho, [], d);
    end
else % Else, start from zero
    rho = zeros(size(rho), 'single');
end

% -------------------------------------------------------------------------
% Select number of channels
idx = randperm(size(x,4));
N = size(x,4);      % < e.g., 2
% N = 8;
x = x(:,:,:,idx(1:N),:);
s = s(:,:,:,idx(1:N),:);

ll  = NaN;
llm = NaN;
llp = NaN;
for it=1:itmax
    
    fprintf('Iteration %d\n', it);
    [rho,llm,llp] = multicoil_sense_mean_map(x, s, A, rho, msk, dir, prm, vs, llm, llp);
    ll = [ll (llm+llp)];
    multicoil_plot_mean(rho, C, ll, vs);
    
    
    % ---------------------------------------------------------------------
    % Check gain
    if it > 1
        gain = (ll(end) - ll(end-1))/(max(ll(2:end), [], 'omitnan') - min(ll(2:end), [], 'omitnan'));
        if verbose > 0
            if gain > 0
                sgn = '+';
            elseif gain < 0
                sgn = '-';
            else
                sgn = '=';
            end
            fprintf('Gain: %20.10g (%s)\n', gain, sgn);
        end
        if abs(gain) < tol
            break
        end
    end
end