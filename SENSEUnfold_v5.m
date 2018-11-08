function Image = SENSEUnfold_v5(FFTImage, estS, UnfInd, accel, precompPinv, weights)
% FORMAT Image = SENSEUnfold_v5(FFTImage, estS, UnfInd, accel, precompPinv, weights)
%
% FFTImage    - Acquired aliased image  (2D/spatial) [nCoil nLine nPart 1] 
% estS        - Estimated sensitivities (2D/spatial) [nCoil nLine nPart 1]
% UnfInd      - Aliased locations (indices/spatial)  [iLine iPart]
% accel       - Acceleration factor [part line]
% precompPinv - g.precomputedPinv [true/false]
% weights     - Weights of aliased locations (spatial) [nPart nLine]
%
% Image       - Unaliased image (spatial) [nLine*NPart 1]


    [~, Nky, Nkz] = size(FFTImage);

    % Unfolding:
    Image = complex(zeros(Nky*Nkz, 1, 'single'));

    l = UnfInd(:,1);
    k = UnfInd(:,2);

    W = diag(inv(diag(weights))); % density + phase correction. Slows the unfolding down a little

    if precompPinv
        nNkz = Nkz/accel(1);
        nNky = Nky/accel(2);
        for zVox = 1: nNkz %Nkz/accel(1)
            for yVox = 1 : nNky %Nky/accel(2)
                samplesYZ = mod(l + yVox - 2, Nky) + (mod(k + zVox - 2, Nkz))*Nky +1;
                Image(samplesYZ) = W.*( estS(:,samplesYZ) ).' * FFTImage(:, yVox, zVox);
            end
        end
    else
        for zVox = 1:Nkz/accel(1)
            for yVox = 1:Nky/accel(2)

                samplesYZ = mod(l + yVox - 2, Nky) + (mod(k + zVox - 2, Nkz))*Nky +1;

                % Unaliase using sensitivity matrix:
                SMatrix = estS(:, samplesYZ);

                % Can use first coil in check since SumOfSquares mask used for all
                vox = SMatrix(1,:) ~= 0;
                if sum(vox) > 0
                    % At least one voxel contributes => do something
                    Image(samplesYZ(vox)) = W(vox).*(SMatrix(:,vox) \ FFTImage(:, yVox, zVox));
                end
            end
        end
    end
end