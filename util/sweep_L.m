function [results,figID] = sweep_L(t,ground_truth,signal,diff_parameters,L_settings,figID)
% This function plots the RMSE (if groun_truth is available) and the
% proposed Differentiator cost for the values of L specified in L_settings.

    % Differentiator parameters
    n_d = diff_parameters{2};
    n_f = diff_parameters{3};
    q = diff_parameters{4};
    trans_diff = diff_parameters{5};

    if isnan(n_d)
        error('Differentiator parameter n_d needs to be selected, NaN is not valid!')
    end

    if isnan(n_f)
        n_f = 0;
    end

    if isnan(trans_diff)
        trans_diff = 0.05*(t(end)-t(1));
    end
    
    % L parameter settings
    exp_min = L_settings(1);
    exp_max = L_settings(2);
    N_L = L_settings(3);

    if N_L < 2
        error('Number of L values to investigate not valid!');
    end

    L_vec = logspace(exp_min,exp_max,N_L);
    
    % Initialisation iteration variables
    results = zeros(3,N_L);
    L = L_vec(1);
    if trans_diff > 0
        if ~isnan(ground_truth)
            [~,indx] = min(abs((t-t(1))-trans_diff));
            ground_truth = ground_truth(indx:end-indx,:);
        end
    end
    
    signal_copy = signal;   % avoids cutting the signal at each iteration

    for i=1:N_L
        diff_output = differentiator(signal_copy,t,n_d,n_f,L,q,0);

        % Discard transient
        if trans_diff > 0
            diff_output = diff_output(indx:end-indx,:);
            signal = signal_copy(indx:end-indx,:);
        end
        
        % Computation of RMSE
        if ~isnan(ground_truth)
            rmse = 0;
            for j=1:length(ground_truth)
                rmse = rmse + (diff_output(j,1)-ground_truth(j))^2;
            end
            results(1,i) = sqrt(rmse/length(ground_truth));
        end

        % Computation of residual statistics
        results(2,i) = mean(diff_output(:,1)-signal);
        results(3,i) = var(diff_output(:,1)-signal);
        
        if i ~= N_L
            L = L_vec(i+1);
        end
        
        % Progression log
        disp(newline); disp(newline); disp(newline);
        disp([num2str(i),'/',num2str(N_L)]);
        disp(newline); disp(newline); disp(newline);
    end
    
    % Plot pointwise error metrics (only if ground truth available)
    if ~isnan(ground_truth)
    
        figID = figID + 1;
        figure(figID);
        min_rmse = min(results(1,:));
        indx_min_rmse = results(1,:)==min_rmse;
        L_min_rmse = L_vec(indx_min_rmse);
        str1 = ['$L=',num2str(round(L_min_rmse(1))),'$'];
        semilogx(L_vec,results(1,:),'LineWidth',2);
        hold on;
        xline(L_min_rmse(1),'-',str1,'LineWidth',2,'interpreter','latex');
        xlabel('$L$','interpreter','latex'); ylabel('RMSE');
        xlim([10^exp_min,10^exp_max]);
        ax = gca;
        ax.FontSize = 25;
        pbaspect([0.8,1,1])
    end

    % Plot proposed Differentiator cost
    figID = figID + 1;
    figure(figID);
    C_diff1 = abs(results(2,:))+sqrt(results(3,:));
    min_cost = min(C_diff1);
    indx_min_cost = C_diff1==min_cost;
    L_min_cost = L_vec(indx_min_cost);
    str2 = ['$L=',num2str(round(L_min_cost(1))),'$'];
    semilogx(L_vec,C_diff1,'LineWidth',2);
    hold on;
    xline(L_min_cost(1),'-',str2,'LineWidth',2,'interpreter','latex');
    xlabel('$L$','interpreter','latex'); ylabel('$\mathcal{C}_{diff}$','interpreter','latex');
    xlim([10^exp_min,10^exp_max]);
    ax = gca;
    ax.FontSize = 25;
    pbaspect([0.8,1,1])

    figID = figID + 1;
    figure(figID);
    C_diff2 = sqrt(results(3,:));
    min_cost = min(C_diff2);
    indx_min_cost = C_diff2==min_cost;
    L_min_cost = L_vec(indx_min_cost);
    str2 = ['$L=',num2str(round(L_min_cost(1))),'$'];
    semilogx(L_vec,C_diff2,'LineWidth',2);
    hold on;
    xline(L_min_cost(1),'-',str2,'LineWidth',2,'interpreter','latex');
    xlabel('$L$','interpreter','latex'); ylabel('$\mathcal{C}_{diff}$','interpreter','latex');
    xlim([10^exp_min,10^exp_max]);
    ax = gca;
    ax.FontSize = 25;
    pbaspect([0.8,1,1])

end

