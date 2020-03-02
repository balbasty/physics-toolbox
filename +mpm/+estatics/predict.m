function pred = predict(out, in, opt)

TEs  = opt.predict(:)';
if isempty(TEs), pred = {}; return; end

pred = in;

cnt = struct;
for i=1:numel(in)
    in1      = in{i};
    
    % ---------------------------------------------------------------------
    % Prepare predicted structure
    nii0           = pred{i}.nii(1);
    pred{i}.nii    =  nifti;
    pred{i}.echoes = {};
    
    % ---------------------------------------------------------------------
    % Read volume info
    ymat = out.mat;                 % Model orientation matrix
    xmat = in1.mat;                 % Observed orientation matrix

    ct  = in1.type;                 % Contrast name
    mat = ymat\xmat;                % Observed-to-Model matrix
    dim = in1.dim;                  % Observed dimensions
    
    if ~isfield(cnt, ct), cnt.(ct) = 0; end
    cnt.(ct) = cnt.(ct) + 1;
    
    % ---------------------------------------------------------------------
    % Load model data
    y0  = utils.pull(single(out.(ct).dat()), mat, dim);
    r2  = utils.pull(single(out.R2s.dat()),  mat, dim);

    
    % ---------------------------------------------------------------------
    % Predict echo
    for j=1:numel(TEs)
        t = TEs(j);
        x = exp(y0 - r2*t);
        pred{i}.nii(j) = nii0;
        fname = sprintf('pred_%s(%d)_TE=%fms.nii', ct, cnt.(ct), t*1E3);
        pred{i}.nii(j).dat.fname = fullfile(opt.out.folder, fname);
        create(pred{i}.nii(j));
        pred{i}.nii(j).dat(:) = x(:);
        pred{i}.echoes{j} = struct('TE', t, 'dat', pred{i}.nii(j).dat, 'var', pred{i}.var);
    end
    
end

