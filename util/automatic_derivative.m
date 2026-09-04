function [t,derivatives_estimate] = automatic_derivative(t,signal,n_derivatives,diff_parameters,SG_parameters,complexity_settings)
% This function automatically estimates the derivatives of signal up to the
% order n_derivatives with the staircase architecture.

    t_init_total = tic;

    derivatives_estimate = zeros(length(signal),n_derivatives+1);

    [diff_parameters,SG_parameters,complexity_settings] = adjust_user_settings(n_derivatives,diff_parameters,SG_parameters,complexity_settings);

    % Extract Differentiator settings
    n_f = diff_parameters{3}(1);
    q = diff_parameters{4}(1);
    trans_diff = diff_parameters{5}(1);

    if isnan(n_f)
        n_f = 0;
    end

    if isnan(trans_diff)
        trans_diff = 0.05*(t(end)-t(1));
    end

    % Extract Savitzky-Golay filter parameters
    m = SG_parameters{1}(1);
    if isnan(m)
        m = 2;
    end

    % Denoise the input signal
    n_d = 0;
    diff_parameters{2}(1) = n_d;
    L = diff_parameters{1}(1);
    if isnan(L)   % L is automatically tuned
        L = autotune_Differentiator(signal,t,diff_parameters,complexity_settings,1);
        diff_parameters{1}(1) = L;
    end

    signal_diff = differentiator(signal,t,n_d,n_f,L,q,0);

    if ~(SG_parameters{2}(1) == 0)    % if window selected is not zero, calculate it
        fl = SG_parameters{2}(1);
        if isnan(fl)    % SG filter window is automatically tuned
            fl = autotune_SG_filter_window(signal,t,diff_parameters,SG_parameters,n_d,complexity_settings,1);
        end
        if mod(fl,2) == 0
            fl = fl + 1;
        end
        signal_denoised = sgolayfilt(signal_diff,m,fl);
    else
        signal_denoised = signal_diff;
    end

    disp('Denoising of signal completed!');
    disp(newline);

    starting_signal = signal_denoised;

    derivatives_estimate(:,1) = starting_signal;

    if n_derivatives > 0
        for i=1:n_derivatives
            n_d = 1;
            diff_parameters{2}(i+1) = n_d;   % estimate only first-derivative at each iteration

            % Estimate first-derivative of previously computed signal
            L = diff_parameters{1}(i+1);
            if isnan(L)
                L = autotune_Differentiator(starting_signal,t,diff_parameters,complexity_settings,i+1);
                diff_parameters{1}(i+1) = L;
            end

            diff_output = differentiator(starting_signal,t,n_d,n_f,L,q,0);    
            
            % Savitzky-Golay filtering of Differentiator first-derivative estimate
            if ~(SG_parameters{2}(i+1) == 0)
                fl = SG_parameters{2}(i+1);
                if isnan(fl)
                    fl = autotune_SG_filter_window(starting_signal,t,diff_parameters,SG_parameters,n_d,complexity_settings,i+1);
                end
                if mod(fl,2) == 0
                    fl = fl + 1;
                end
                derivative = sgolayfilt(diff_output(:,2),m,fl);

            else

                derivative = diff_output(:,2);
            end
            
            derivatives_estimate(:,i+1) = derivative;
            if i < n_derivatives
                starting_signal = derivative;
            end
    
            disp(['Derivative of order ',num2str(i),' estimated with staircase method!']);
            disp(newline);
        end

    end

    % Elimination differentiator transitory
    if trans_diff > 0
        [~,indx] = min(abs((t-t(1))-trans_diff));
        derivatives_estimate = derivatives_estimate(indx:end-indx,:);
        t = t(indx:end-indx);
    end

    t_fin_total = toc(t_init_total);
    disp(['Total time to estimate all the derivatives: ',num2str(t_fin_total),' [s]']);

end