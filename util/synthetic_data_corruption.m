function corrupted_signal = synthetic_data_corruption(signal_to_corrupt,noise_type,noise_var,additive_noise)
% This function synthetically adds noise of the type noise_type, in the way
% specified by "additiv_noise", with characteristic specified by noise_var,
% to the signal signal_to_corrupt.

    switch noise_type
        case 0
            eta = zeros(length(signal_to_corrupt),1);
        case 1
            rng(3,'twister')   % for repeatability
            eta = wgn(length(signal_to_corrupt),1,noise_var,'linear');
        otherwise
            error('Noise type not valid.');
    end
    
    % Measured time-series data
    switch additive_noise
        case 0
            corrupted_signal = signal_to_corrupt.*(eta+ones(size(eta,1),size(eta,2)));
        case 1
            corrupted_signal = signal_to_corrupt + eta;
        otherwise
            error('Data corruption mode not valid.');
    end

end