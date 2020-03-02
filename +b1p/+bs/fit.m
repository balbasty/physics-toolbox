function model = fit(dat,opt)


[sens_pos, meanim_pos] = b1m.fit_robust(dat.coils(:,:,:,:,1), 'VoxelSize', dat.vs, 'Verbose', 2);
[sens_neg, meanim_neg] = b1m.fit_robust(dat.coils(:,:,:,:,2), 'VoxelSize', dat.vs, 'Verbose', 2);



% --- Define option and model structures
if nargin < 2, opt = struct; end
opt   = b1p.bs.opt(opt);
dim   = [size(dat.coils) 1];
dim   = dim(1:3);
s      = struct;
s.mean = zeros(dim, 'like', single(1i));
s.sens = zeros([dim size(dat.coils, 4)], 'like', single(1i));
s.b1   = zeros(dim, 'single');
s.vs   = vs;
model.sens(:,:,:,:) = 1;
if ~opt.b1.log, model.b1(:,:,:) = 1; end

% --- Initialise model
for it=1:5
    model = b1p.bs.mean.update(dat,model,opt);
    model = b1p.bs.sens.init(dat,model,opt);
    model = b1p.bs.b1.init(dat,model,opt);
end

% --- Variational optimisation
ll = [];
lls1 = 0;
llm1 = 0;
llb1 = 0;
for it=1:opt.nbiter
    
    [model,llx1,lls1] = b1p.bs.sens.update(dat,model,opt);
    ll = [ll (llx1 + lls1 + llm1 + llb1)];
    plot(ll); drawnow
    [model,llx1,llm1] = b1p.bs.mean.update(dat,model,opt);
    ll = [ll (llx1 + lls1 + llm1 + llb1)];
    plot(ll); drawnow
    [model,llx1,llb1] = b1p.bs.b1.update(dat,model,opt);
    ll = [ll (llx1 + lls1 + llm1 + llb1)];
    plot(ll); drawnow
    
end