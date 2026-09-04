function fl_opt = autotune_SG_filter_window(signal,t,diff_parameters,SG_parameters,n_derivative,complexity_settings,iteration)
% This function automatically computes the SG filter window size based on
% the SG filter cost and the persistence of the frequency maximizing the
% PSD of the residual criterion.

    % Auxiliary parameters
    t_init = tic;

    use_SG_cost = complexity_settings{6}(iteration);

    % Savitzky-Golay filter settings
    m = SG_parameters{1}(iteration);
    SG_n_points_min = SG_parameters{3}(iteration);
    SG_n_points_max = SG_parameters{4}(iteration);
    delta_points = SG_parameters{5}(iteration);
    lambda = SG_parameters{6}(iteration);
    n_windows = SG_parameters{7}(iteration);
    overlap = SG_parameters{8}(iteration);
    epsilon = SG_parameters{9}(iteration);

    if isnan(epsilon)
        epsilon = 0.05;
    end

    % Differentiator parameters
    L_opt = diff_parameters{1}(iteration);
    n_d = diff_parameters{2}(iteration);
    n_f = diff_parameters{3}(iteration);
    q = diff_parameters{4}(iteration);
    trans_diff = diff_parameters{5}(iteration);

    if isnan(n_f)
        n_f = 0;
    end

    if isnan(trans_diff)
        trans_diff = 0.05*(t(end)-t(1));
    end

    % Compute the representative chattering signal
    L_max = (1+epsilon)*L_opt;

    diff_output_opt = differentiator(signal,t,n_d,n_f,L_opt,q,0);
    diff_output_max = differentiator(signal,t,n_d,n_f,L_max,q,0);

    if trans_diff > 0
        [~,indx] = min(abs((t-t(1))-trans_diff));
        diff_output_opt = diff_output_opt(indx:end,:);
        diff_output_max = diff_output_max(indx:end,:);
    end

    chattering = diff_output_max(:,n_derivative+1) - diff_output_opt(:,n_derivative+1);

    diff_output_opt = diff_output_opt(:,n_derivative+1); % signal of interest

    % Set to default the parameters not specified by the user
    if isnan(m)
        m = 2;
    end

    if isnan(lambda)
        lambda = 0.5;
    end

    if isnan(n_windows)
        n_windows = [2,4,8,16];
    end

    if isnan(overlap)
        overlap = 0.5*ones(1,length(n_windows));
    end

    if isnan(delta_points)
        if use_SG_cost
            delta_points = 6;
        else
            delta_points = 60;
        end
        
    end

    if isnan(SG_n_points_min)
        SG_n_points_min = m+1;
    end

    if mod(SG_n_points_min,2) == 0
        SG_n_points_min = SG_n_points_min + 1;
    end

    if isnan(SG_n_points_max)
        SG_n_points_max = 20001;
        if length(chattering)<SG_n_points_max
            SG_n_points_max = length(chattering)-1;
        end
    end

    if mod(SG_n_points_max,2) == 0
        SG_n_points_max = SG_n_points_max + 1;
    end

    %% Find window that minimises the SG filter cost
    if use_SG_cost
        smoothed_chat_min = sgolayfilt(chattering,m,SG_n_points_min);
        norm_factor_chat = 1/(length(smoothed_chat_min)-1)*sum(abs(diff(smoothed_chat_min)));
        smoothed_diff_max = sgolayfilt(diff_output_opt,m,SG_n_points_max);
        norm_factor_dist = 1/var(smoothed_diff_max);
        % norm_factor_dist = 1/sqrt(var(smoothed_diff_max));
    
        cost = @(x) cost_window_SG_filter(x,m,lambda,t,chattering,diff_output_opt,norm_factor_chat,norm_factor_dist);
    
        % Find optimal L value that minimises the residual variance with simulated annealing
        max_iterations_SG = complexity_settings{3}(iteration);
        if isnan(max_iterations_SG)
            max_iterations_SG = 50;
        end
        tolerance_SG = complexity_settings{4}(iteration);
        if isnan(tolerance_SG)
            tolerance_SG = 1e-4;
        end
        SA_options = optimoptions('simulannealbnd','MaxIterations',max_iterations_SG,'FunctionTolerance',tolerance_SG,'Display','off');
        if complexity_settings{7}(iteration)
         rng(123,'twister')   % for repeatability
        end
        window_min_cost = simulannealbnd(cost,t(SG_n_points_min)-t(1),t(SG_n_points_min)-t(1),t(SG_n_points_max)-t(1),SA_options);
    
        disp(['Window for Savitzky-Golay filter that minimises the cost found is ',num2str(round(window_min_cost,2,'significant'))]);
    
        [~,window_min_cost] = min(abs((t-t(1))-window_min_cost));
        if mod(window_min_cost,2) == 0
            window_min_cost = window_min_cost + 1;
        end
    else
        window_min_cost = SG_n_points_max;
    end

    %% Remove interval of persistent max PSD frequency at cost minimiser
    use_persistence = complexity_settings{5}(iteration);
    
    if use_persistence
        smoothed_diff = sgolayfilt(diff_output_opt,m,window_min_cost);
    
        freq_at_cost_min = zeros(1,length(n_windows));
        dt=t(2)-t(1);
        parfor i=1:length(n_windows)
            [~,freqVec,~,psd] = spectrogram(smoothed_diff-diff_output_opt,round(length(smoothed_diff)/n_windows(i)),round(overlap(i)*length(smoothed_diff)/n_windows(i)),[],1/dt);
            meanPSD = mean(psd,2);
            tmp = db(meanPSD,"power");
            [~,indx_psd] = max(tmp);
            freq_at_cost_min(i) = freqVec(indx_psd);
        end
    
        persistence = true(1,length(n_windows));
        fl_opt_vec = zeros(1,length(n_windows));
        fl = window_min_cost;
    
        while any(persistence)
            fl = fl-delta_points;

            if fl <= m
                fl_opt = m+1;
                break;
            end
    
            if mod(fl,2) == 0
                fl = fl - 1;
            end
            % Apply SG filter to the chattering component and optimal Differentiator output
            smoothed_diff = sgolayfilt(diff_output_opt,m,fl);
    
            parfor i=1:length(n_windows)
                % Compute PSD of the residue and extract frequency of max power
                [~,freqVec,~,psd] = spectrogram(smoothed_diff-diff_output_opt,round(length(smoothed_diff)/n_windows(i)),round(overlap(i)*length(smoothed_diff)/n_windows(i)),[],1/dt);
                meanPSD = mean(psd,2);
                tmp = db(meanPSD,"power");
                [~,indx_psd] = max(tmp);
                max_psd = freqVec(indx_psd);
        
                if max_psd ~= freq_at_cost_min(i)
                    persistence(i) = false;
                    fl_opt_vec(i) = fl;
                end
            end
        end
    
        if ~exist('fl_opt','var')
            fl_opt = round(mean(fl_opt_vec));
        end
    else
        fl_opt = window_min_cost;
    end

    if mod(fl_opt,2) == 0
        fl_opt = fl_opt + 1;
    end

    if fl_opt < SG_n_points_min
        fl_opt = SG_n_points_min;
    end

    t_fin = toc(t_init);
    disp(['Savitzky-Golay filter window calculation: ',num2str(round(t_fin,3,'significant')),' [s]']);
    if use_persistence
        disp(['Estimated final window is ',num2str(round(t(fl_opt)-t(1),2,'significant'))]);
    else
        disp('PSD information not used, final window is the window that minimises the cost');
    end

end