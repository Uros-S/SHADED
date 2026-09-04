function [results_rmse,SG_cost,freq_max_psd,figID] = sweep_SG_window(t,ground_truth,signal,fl_settings,diff_parameters,SG_settings,figID)
% This function plots the normalized RMSE (if groun_truth is available) and 
% the proposed SG filter cost for the values of window sizes specified in 
% fl_settings and with the selected Differentiator gain L in 
% diff_parameters. It also plots the frequencies that maximize the PSD of
% the residual w.r.t. the SG filter window size length for every number of
% windows selected for the Welch method.
    
    dt = t(2)-t(1);     % assuming constant sampling rate
    
    % Differentiator parameters
    L = diff_parameters{1};
    n_d = diff_parameters{2};
    n_f = diff_parameters{3};
    q = diff_parameters{4};
    trans_diff = diff_parameters{5};

    if isnan(L)
        error('Differentiator gain L needs to be selected, NaN is not valid!')
    end

    if isnan(n_d)
        error('Differentiator parameter n_d needs to be selected, NaN is not valid!')
    end

    if isnan(n_f)
        n_f = 0;
    end

    if isnan(trans_diff)
        trans_diff = 0.05*(t(end)-t(1));
    end
    
    % Savitzky-Golay filter settings
    m = SG_settings{1};    
    lambda = SG_settings{6}(1);
    n_windows = SG_settings{7}(1);
    overlap = SG_settings{8}(1);
    epsilon = SG_settings{9}(1);
    length_min = fl_settings(1);
    length_max = fl_settings(2);
    N_SG = fl_settings(3);

    if isnan(m)
        m = 2;
        %error('SG filter order m needs to be selected, NaN is not valid!')
    end

    if isnan(epsilon)
        epsilon = 0.05;
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

    if  length(t) <= 20001
        SG_settings{4}(1) = length(t);
        if mod(SG_settings{4}(1),2) == 0    % SG filter window size must be odd
            SG_settings{4}(1) = SG_settings{4}(1) - 1;
        end
    end

    if N_SG < 2
        error('Number of SG filter window lengths values to investigate not valid!');
    end

    if length_min <= m || length_max >= length(t) || length_min >= length_max
        error('Maximum and minimum window lengths values to investigate not valid!');
    end

    if isnan(length_min)
        length_min = m+1;
        fl_settings(1) = length_min;
    end

    if isnan(length_max)
        length_max = 20001;
        fl_settings(2) = length_max;
    end

    length_min = length_min*dt;
    length_max = length_max*dt;

    SG_lengths_vec_original = linspace(length_min, length_max, N_SG);
    SG_lengths_vec = zeros(1,N_SG);
    for i=1:N_SG    % convert seconds to vector index
        SG_lengths_vec(i) = round(SG_lengths_vec_original(i)/dt);
        if mod(SG_lengths_vec(i),2) == 0    % SG filter window size must be odd
            SG_lengths_vec(i) = SG_lengths_vec(i) + 1;
        end
    end
    
    % Initialisation iteration variables
    results_rmse = zeros(1,N_SG);
    TV_chat = zeros(1,N_SG);
    D_sig = zeros(1,N_SG);
    SG_cost = zeros(1,N_SG);
    freq_max_psd = zeros(length(n_windows),N_SG);

    fl = SG_lengths_vec(1);
    if trans_diff > 0
        if ~isnan(ground_truth)
            ground_truth = ground_truth(round(trans_diff/dt):end,:);
        end
    end

    % Compute representative chattering signal
    L_max = (1+epsilon)*L;

    diff_output_opt = differentiator(signal,t,n_d,n_f,L,q,0);
    diff_output_max = differentiator(signal,t,n_d,n_f,L_max,q,0);

    if trans_diff > 0
        [~,indx] = min(abs((t-t(1))-trans_diff));
        diff_output_opt = diff_output_opt(indx:end,:);
        diff_output_max = diff_output_max(indx:end,:);
    end

    chattering = diff_output_max(:,1) - diff_output_opt(:,1);

    % Compute other quantities needed for the SG cost
    if mod(fl_settings(1),2) == 0
        fl_settings(1) = fl_settings(1) + 1;
    end

    if mod(fl_settings(2),2) == 0
        fl_settings(2) = fl_settings(2) + 1;
    end

    smoothed_chat_min = sgolayfilt(chattering,m,fl_settings(1));
    norm_factor_chat = 1/(length(smoothed_chat_min)-1)*sum(abs(diff(smoothed_chat_min)));
    if ~isnan(SG_settings{4}(1))
        smoothed_diff_max = sgolayfilt(diff_output_opt(:,1),m,SG_settings{4}(1));
    else
        smoothed_diff_max = sgolayfilt(diff_output_opt(:,1),m,20001);
    end
    norm_factor_dist = 1/var(smoothed_diff_max);

    %% Iterate over the selected windows
    for i=1:N_SG
        diff_output = differentiator(signal,t,n_d,n_f,L,q,0);
        diff_output_filt = sgolayfilt(diff_output(:,1),m,fl);

        % Discard transient
        if trans_diff > 0
            diff_output_filt = diff_output_filt(round(trans_diff/dt):end,:);
        end
        
        % Computation of RMSE
        if ~isnan(ground_truth)
            rmse = 0;           % Root Mean Squared Error
            for j=1:length(ground_truth)
                rmse = rmse + (diff_output_filt(j,1)-ground_truth(j))^2;
            end
    
            results_rmse(1,i) = sqrt(rmse/length(ground_truth));
        end

        % Compute SG cost
        smoothed_chat = sgolayfilt(chattering,m,fl);
        smoothed_diff = sgolayfilt(diff_output_opt(:,1),m,fl);

        cost_chat = 1/(length(smoothed_chat)-1)*sum(abs(diff(smoothed_chat)));
        cost_dist = 1/var(smoothed_diff);

        cost_chat = 1/norm_factor_chat*cost_chat;
        cost_dist = 1/norm_factor_dist*cost_dist;

        TV_chat(i) = cost_chat;
        D_sig(i) = cost_dist;
        SG_cost(i) = (1-lambda)*cost_chat + lambda*cost_dist;

        % Compute frequency at the PSD maximum
        for k=1:length(n_windows)
            [~,freqVec,~,psd] = spectrogram(smoothed_diff(:,1)-diff_output_opt(:,1),round(length(smoothed_diff)/n_windows(k)),round(overlap(k)*length(smoothed_diff)/n_windows(k)),[],1/dt);
            meanPSD = mean(psd,2);
            tmp = db(meanPSD,"power");
            [~,indx_psd] = max(tmp);
    
            freq_max_psd(k,i) = freqVec(indx_psd);
        end
        
        if i ~= N_SG
            fl = SG_lengths_vec(i+1);
        end
        
        % Progression log
        disp(newline); disp(newline); disp(newline);
        disp([num2str(i),'/',num2str(N_SG)]);
        disp(newline); disp(newline); disp(newline);
    end
    
    %% Plot RMSE
    if ~isnan(ground_truth)
    
        figID = figID + 1;
        figure(figID);
        min_rmse = min(results_rmse);
        indx_min_rmse = results_rmse==min_rmse;
        norm_factor_rmse = 0;
        diff_output = differentiator(signal,t,n_d,n_f,L,q,0);
        if trans_diff > 0
            diff_output = diff_output(round(trans_diff/dt):end,:);
        end
        for j=1:length(ground_truth)
            norm_factor_rmse = norm_factor_rmse + (diff_output(j,1)-ground_truth(j))^2;
        end
        norm_factor_rmse = sqrt(norm_factor_rmse/length(ground_truth));

        SG_length_min_rmse = SG_lengths_vec_original(indx_min_rmse);
        str1 = ['$w^*=',num2str(round(SG_length_min_rmse(1),3,'significant')),'$'];
        plot(SG_lengths_vec_original,1/norm_factor_rmse*results_rmse(1,:),'LineWidth',2);
        hold on;
        xline(SG_length_min_rmse(1),'-',str1,'LineWidth',2,'interpreter','latex');
        xlabel('SG window length [s]'); 
        ylabel('Normalized RMSE');
        title(['$L=',num2str(L),'$'],'interpreter','latex');
        xlim([length_min,length_max]);
        ax = gca;
        ax.FontSize = 25;
        pbaspect([0.8,1,1])

    end

    %% Plot SG cost and its components
    figID = figID + 1;
    figure(figID);
    plot(SG_lengths_vec_original,TV_chat,'LineWidth',2);
    xlabel('SG window length [s]');
    ylabel('$TV_{chat}$','interpreter','latex');
    title(['$L=',num2str(L),'$'],'interpreter','latex');
    xlim([length_min,length_max]);
    ax = gca;
    ax.FontSize = 25;
    pbaspect([0.8,1,1])

    figID = figID + 1;
    figure(figID);
    plot(SG_lengths_vec_original,D_sig,'LineWidth',2);
    xlabel('SG window length [s]');
    ylabel('$D_{sig}$','interpreter','latex');
    title(['$L=',num2str(L),'$'],'interpreter','latex');
    xlim([length_min,length_max]);
    ax = gca;
    ax.FontSize = 25;
    pbaspect([0.8,1,1])

    figID = figID + 1;
    figure(figID);
    min_SG_cost = min(SG_cost);
    indx_min_SG_cost = SG_cost==min_SG_cost;
    SG_length_min_SG_cost = SG_lengths_vec_original(indx_min_SG_cost);
    str2 = ['$w=',num2str(round(SG_length_min_SG_cost(1),3,'significant')),'$'];
    plot(SG_lengths_vec_original,SG_cost,'LineWidth',2);
    hold on;
    xline(SG_length_min_SG_cost,'-',str2,'LineWidth',2,'interpreter','latex');
    xlabel('SG window length [s]');
    ylabel('$C_{SG}$','interpreter','latex');
    title(['$L=',num2str(L),'$'],'interpreter','latex');
    xlim([length_min,length_max]);
    ax = gca;
    ax.FontSize = 25;
    pbaspect([0.8,1,1])

    %% Plot frequencies that maximize PSD
    for k=1:length(n_windows)
        figID = figID + 1;
        figure(figID);
        plot(SG_lengths_vec_original,freq_max_psd(k,:),'LineWidth',2);
        hold on;
        xline(SG_length_min_SG_cost,'-','LineWidth',2,'interpreter','latex');
        xlabel('SG window length [s]');
        ylabel('\omega_{\text{max}}','interpreter','latex');
        title(['$L=',num2str(L),'$'],'interpreter','latex');
        ax = gca;
        ax.FontSize = 25;
        set(gca,'Yscale','log')
    end

end