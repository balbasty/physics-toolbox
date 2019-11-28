function [sig,mu] = estimate_snr(dat)

    N = size(dat,4);
    ndigits = ceil(log10(N))+1;
    sig  = zeros(1,N);
    mu = zeros(1,N);

    fprintf(['Estimate SNR: %' num2str(ndigits) 's'], '');
    for n=1:N

        fprintf([repmat('\b', [1 ndigits]) '%' num2str(ndigits) 'd'], n);

        % -----------------------------------------------------------------
        % Load one coil image
        xn = single(dat(:,:,:,n));
        xn = reshape(xn, [], 1);

        % -----------------------------------------------------------------
        % Compute histogram
        xmax = max(xn,[],'omitnan');
        xmin = min(xn,[],'omitnan');
        cn = linspace(xmin,xmax,129);
        cn = (cn(2:end)+cn(1:end-1))/2;
        bwn = (xmax - xmin)/128;
        [xn,wn] = spm_imbasics('hist', double(xn), cn(:), 'KeepZero', false);
        clear c1

        % -----------------------------------------------------------------
        % Fit Rician mixture
        [PI,NU,SIG] = spm_rice_mixture(double(wn), double(xn), 2);
        if NU(1) == 0 && NU(2) == 0
            % If failed, fit Gaussian mixture
            [~,MU,A,PI] = spm_gmm(double(xn), 2, double(wn),...
                'GaussPrior', {[],10,[],10},'BinWidth',bwn);
            if MU(1) <= MU(2)
                sig(n)  = 1./A(1);
                mu(n) = MU(2);
            else
                sig(n)  = 1./A(2);
                mu(n) = MU(1);
            end
        else
            if NU(1) <= NU(2)
                sig(n)  = SIG(1)^2;
                mu(n) = NU(2);
            else
                sig(n)  = SIG(2)^2;
                mu(n) = NU(1);
            end
        end

    end
    fprintf('\n');

    sig = exp(mean(log(sig)));
    
end