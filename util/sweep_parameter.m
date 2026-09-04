function [results_param,figID] = sweep_parameter(t,signal,diff_parameters,SG_parameters,complexity_settings,param_settings,param_number,figID)

    dt = t(2)-t(1);

    param_min = param_settings(1);
    param_max = param_settings(2);
    N_param = param_settings(3);

    if param_number == 3
        param_vec = round(linspace(param_min, param_max, N_param));
    else
        param_vec = linspace(param_min, param_max, N_param);
    end
    results_param = zeros(2,N_param);

    n_derivative = 0;   % consider denoising case

    for i=1:N_param
        if param_number == 1
            SG_parameters{6} = param_vec(i);
        elseif param_number == 2
            SG_parameters{9} = param_vec(i);
        else
            if mod(param_vec(i),2) == 0
                param_vec(i) = param_vec(i) + 1;
            end
            SG_parameters{4} = param_vec(i);
        end

        complexity_settings{5} = false;      % Consider w'
        results_param(1,i) = autotune_SG_filter_window(signal,t,diff_parameters,SG_parameters,n_derivative,complexity_settings,1);

        complexity_settings{5} = true;       % Consider \hat{w}
        results_param(2,i) = autotune_SG_filter_window(signal,t,diff_parameters,SG_parameters,n_derivative,complexity_settings,1);

        % Progression log
        disp(newline); disp(newline); disp(newline);
        disp([num2str(i),'/',num2str(N_param)]);
        disp(newline); disp(newline); disp(newline);
    end

    figID = figID + 1;
    figure(figID);
    plot(param_vec,dt*results_param(1,:),'LineWidth',2);
    if param_number == 1
        xlabel('$\lambda$','interpreter','latex');
    elseif param_number == 2
        xlabel('$\epsilon$','interpreter','latex');
    else
        xlabel('$w_{max}$','interpreter','latex');
    end
    ylabel('$w$ prime [s]','interpreter','latex');
    title(['$L=',num2str(diff_parameters{1}),'$'],'interpreter','latex');
    ax = gca;
    ax.FontSize = 25;

    figID = figID + 1;
    figure(figID);
    plot(param_vec,dt*results_param(2,:),'LineWidth',2);
    if param_number == 1
        xlabel('$\lambda$','interpreter','latex');
    elseif param_number == 2
        xlabel('$\epsilon$','interpreter','latex');
    else
        xlabel('$w_{max}$','interpreter','latex');
    end
    ylabel('$w$ hat [s]','interpreter','latex');
    title(['$L=',num2str(diff_parameters{1}),'$'],'interpreter','latex');
    ax = gca;
    ax.FontSize = 25;

end