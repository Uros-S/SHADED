function figID = plot_signals_theoretical_model(model,t,s,signal,diff_output,diff_output_guess,model_param,tau,figID)
% This function the results for the theoretical models, comparing the
% automated procedure versus the method with Guessed parameters

    if model == 1
        sigma = model_param(1);
        rho = model_param(2);
        beta = model_param(3);
    elseif model == 2
        x = s;
        a = model_param(1);
        b = model_param(2);
        c = model_param(3);
        d = model_param(4);
        r = model_param(5);
        s = model_param(6);
        I = model_param(7);
        x_R = model_param(8);
    elseif model == 3
        a = model_param(1);
        b = model_param(2);
        I_ext = model_param(3);
        tau = model_param(4);
    elseif model == 4
        C      = model_param(1);
        A      = model_param(2);
        B      = model_param(3);
        a      = model_param(4);
        b      = model_param(5);
        V_0    = model_param(6); 
        nu_max = model_param(7);
        r      = model_param(8); 
        p      = model_param(9);
    end

    %% Plot signals
    % Measured and noise-free signal time-series
    figID = figID + 1;
    figure(figID);
    plot(t,signal,'LineWidth',1,'Color',[0,0.4470,0.7410,0.7]);
    hold on;
    if model == 1 || model == 3 || model == 5
        plot(t,s(:,1),'k','LineWidth',1.5);
    elseif model == 2
        plot(t,x(:,1),'k','LineWidth',1.5);
    elseif model == 4
        plot(t,s(:,2)-s(:,3),'k','LineWidth',1.5);
    end
    legend('Measured signal','Noise-free signal');
    xlim([t(1),t(end)]);
    pbaspect([2,1,1])
    ax = gca;
    ax.FontSize = 35;

    % Zeroth-order derivative estimates
    figID = figID + 1;
    figure(figID);
    plot(t,diff_output(:,1),'LineWidth',2);
    hold on;
    plot(t,diff_output_guess(:,1),'LineWidth',2);
    hold on;
    if model == 1 || model == 3 || model == 5
        zero_der = s(:,1);
        plt = plot(t,zero_der,'k','LineWidth',2);
    elseif model == 2
        zero_der = x(:,1);
        plt = plot(t,zero_der,'k','LineWidth',2);
    elseif model == 4
        zero_der = s(:,2)-s(:,3);
        plt = plot(t,zero_der,'k','LineWidth',2);
    end
    lgd = legend('Automatic','Guessing','Ground truth');
    set(lgd, 'interpreter', 'latex','LineWidth',2)
    uistack(plt,'top');
    xlim([t(1),t(end)]);
    ax = gca;
    ax.FontSize = 35;
    title('Zeroth-derivative estimate');

    % First-order derivative estimate
    if size(diff_output,2) >= 2
        figID = figID + 1;
        figure(figID);
        plot(t,diff_output(:,2),'LineWidth',2);
        hold on;
        plot(t,diff_output_guess(:,2),'LineWidth',2);
        hold on;
        if model == 1
            first_der = -sigma*(s(:,1)-s(:,2));
            plt = plot(t,first_der,'k','LineWidth',2);
        elseif model == 2
            first_der = -a*x(:,1).^3+b*x(:,1).^2+x(:,2)-x(:,3)+I*ones(length(x(:,1)),1);
            plt = plot(t,first_der,'k','LineWidth',2);
        elseif model == 3
            first_der = s(:,1)-s(:,1).^3/3-s(:,2)+ones(length(s(:,1)),1)*I_ext;
            plt = plot(t,first_der,'k','LineWidth',2);
        elseif model == 4
            first_der = s(:,5)-s(:,6);
            plt = plot(t,first_der,'k','LineWidth',2);
        end
        lgd = legend('Automatic','Guessing','Ground truth');
        set(lgd, 'interpreter', 'latex','LineWidth',2)
        title('First-derivative estimate');
        uistack(plt,'top');
        xlim([t(1),t(end)]);
        ax = gca;
        ax.FontSize = 35;
    end

    % Secon-order derivative estimate
    if size(diff_output,2) >= 3
        figID = figID + 1;
        figure(figID);
        plot(t,diff_output(:,3),'LineWidth',2);
        hold on;
        plot(t,diff_output_guess(:,3),'LineWidth',2);
        hold on;
        if model == 1
            second_der = (s(:,1)-s(:,2))*sigma^2 + (s(:,1).*(ones(length(s(:,1)),1)*rho-s(:,3))-s(:,2))*sigma;
            plt = plot(t,second_der,'k','LineWidth',2);
        elseif model == 2
            second_der = ones(length(x(:,1)),1)*c-x(:,2)+r*(x(:,3)-s*(x(:,1)-ones(length(x(:,1)),1)*x_R))-d*x(:,1).^2+(-3*a*x(:,1).^2+2*b*x(:,1)).*(-a*x(:,1).^3+b*x(:,1).^2+x(:,2)-x(:,3)+ones(length(x(:,1)),1)*I);
            plt = plot(t,second_der,'k','LineWidth',2);
        elseif model == 3
            second_der = -(ones(length(s(:,1)),1)*a+s(:,1)-b*s(:,2))/tau-(s(:,1).^2-ones(length(s(:,1)),1)).*(-s(:,1).^3/3+s(:,1)+ones(length(s(:,1)),1)*I_ext-s(:,2));
            plt = plot(t,second_der,'k','LineWidth',2);
        elseif model == 4
            second_der = zeros(length(s(:,1)),1);
            for i=1:length(s(:,1))
                second_der(i) = 2*b*s(i,6) - 2*a*s(i,5) - a^2*s(i,2) + b^2*s(i,3) + A*a*(p + (4*C*nu_max)/(5*(exp(r*(V_0 - C*s(i,1))) + 1))) - (B*C*b*nu_max)/(4*(exp(r*(V_0 - (C*s(i,1))/4)) + 1));
            end
            plt = plot(t,second_der,'k','LineWidth',2);
        end
        lgd = legend('Automatic','Guessing','Ground truth');
        set(lgd, 'interpreter', 'latex','LineWidth',2)
        title('Second-derivative estimate')
        uistack(plt,'top');
        xlim([t(1),t(end)]);
        ax = gca;
        ax.FontSize = 35;
    end

    % Third-order derivative estimate
    if size(diff_output,2) >= 4
        figID = figID + 1;
        figure(figID);
        plot(t,diff_output(:,4),'LineWidth',2);
        hold on;
        plot(t,diff_output_guess(:,4),'LineWidth',2);
        hold on
        if model == 1
            plt = plot(t,(s(:,2)-s(:,1).*(ones(length(s(:,1)),1)*rho-s(:,3)))*(sigma^2+sigma)-sigma*(sigma^2+(ones(length(s(:,1)),1)*rho-s(:,3))*sigma).*(s(:,1)-s(:,2))+sigma*s(:,1).*(beta*s(:,3)-s(:,1).*s(:,2)),'k','LineWidth',2);
        elseif model == 2
            plt = plot(t,(3*a*x(:,1).^2-2*b*x(:,1)+ones(length(x(:,1)),1)).*(d*x(:,1).^2-ones(length(x(:,1)),1)*c+x(:,2))-(2*d*x(:,1)+ones(length(x(:,1)),1)*r*s-(-3*a*x(:,1).^2+2*b*x(:,1)).^2-(ones(length(x(:,1)),1)*2*b-6*a*x(:,1)).*(-a*x(:,1).^3+b*x(:,1).^2+x(:,2)-x(:,3)+ones(length(x(:,1)),1)*I)).*(-a*x(:,1).^3+b*x(:,1).^2+x(:,2)-x(:,3)+ones(length(x(:,1)),1)*I)-r*(x(:,3)-s*(x(:,1)-ones(length(x(:,1)),1)*x_R)).*(3*a*x(:,1).^2-2*b*x(:,1)+ones(length(x(:,1)),1)*r),'k','LineWidth',2);
        elseif model == 3

        elseif model == 4

        end
        lgd = legend('Diff + SG filter','Competing methodology','Ground truth');
        set(lgd, 'interpreter', 'latex','LineWidth',2)
        title('Third-derivative estimate')
        uistack(plt,'top');
        xlim([t(1),t(end)]);
        ax = gca;
        ax.FontSize = 35;
    end

    % Fourth-order derivative estimate
    if size(diff_output,2) >= 5
        figID = figID + 1;
        figure(figID);
        plot(t,diff_output(:,5),'LineWidth',2);
        hold on;
        plot(t,diff_output_guess(:,5),'LineWidth',2);
        hold on;
        if model == 1
            plt = plot(t,sigma^2*(s(:,1)-s(:,2)).*(ones(length(s(:,1)),1)*rho-s(:,3)-beta*s(:,3)+ones(length(s(:,1)),1)*2*rho*sigma-2*sigma*s(:,3)+2*s(:,1).*s(:,2)+sigma^2)-sigma*(beta*s(:,3)-s(:,1).*s(:,2)).*(s(:,1)+beta*s(:,1)+2*sigma*s(:,1)-sigma*s(:,2))-sigma*(s(:,2)-rho*s(:,1)+s(:,1).*s(:,3)).*(ones(length(s(:,1)),1)*(sigma + rho*sigma+sigma^2+ 1)-sigma*s(:,3)-s(:,1).^2),'k','LineWidth',2);
        elseif model == 2

        elseif model == 3

        elseif model == 4

        end
        lgd = legend('Diff + SG filter','Competing methodology','Ground truth');
        set(lgd, 'interpreter', 'latex','LineWidth',2)
        title('Fourth-derivative estimate')
        uistack(plt,'top');
        xlim([t(1),t(end)]);
        ax = gca;
        ax.FontSize = 35;
    end

    %% Plot differential embeddings
    if size(diff_output,2) >= 2
        figID = figID + 1;
        figure(figID);
        subplot(1,3,1); 
        scatter(zero_der,first_der,0.5,[0,0,0]);
        xlabel('$x_1$','interpreter','latex'); ylabel('$\dot{x}_1$','interpreter','latex');
        title(['Noise-free',newline,'dynamics']);
        ax = gca;
        ax.FontSize = 25;

        subplot(1,3,2); scatter(diff_output_guess(:,1),diff_output_guess(:,2),0.7,[0,0.4470,0.7410]);
        xlabel('$x_{1,D}$','interpreter','latex'); ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        title(['Differentiator reconstruction',newline,'with guess parameters']);
        ax = gca;
        ax.FontSize = 25;

        subplot(1,3,3); scatter(diff_output(:,1),diff_output(:,2),0.7,[0,0.4470,0.7410]);
        xlabel('$x_{1,D}$','interpreter','latex'); ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        title(['Differentiator reconstruction',newline,'automatic']);
        ax = gca;
        ax.FontSize = 25;
    end

    if size(diff_output,2) >= 3
        figID = figID + 1;
        figure(figID);
        subplot(1,3,1); 
        scatter3(zero_der,first_der,second_der,0.5,[0,0,0]);
        xlabel('$x_1$','interpreter','latex'); ylabel('$\dot{x}_1$','interpreter','latex'); zlabel('$\ddot{x}_1$','interpreter','latex');
        title(['Noise-free',newline,'dynamics']);
        ax = gca;
        ax.FontSize = 25;

        subplot(1,3,2); 
        scatter3(diff_output_guess(:,1),diff_output_guess(:,2),diff_output_guess(:,3),0.7,[0,0.4470,0.7410]);
        xlabel('$x_{1,D}$','interpreter','latex'); ylabel('$\dot{x}_{1,D}$','interpreter','latex'); zlabel('$\ddot{x}_{1,D}$','interpreter','latex');
        title(['Differentiator reconstruction',newline,'with guess parameters']);
        ax = gca;
        ax.FontSize = 25;

        subplot(1,3,3); 
        scatter3(diff_output(:,1),diff_output(:,2),diff_output(:,3),0.7,[0,0.4470,0.7410]);
        xlabel('$x_{1,D}$','interpreter','latex'); ylabel('$\dot{x}_{1,D}$','interpreter','latex'); zlabel('$\ddot{x}_{1,D}$','interpreter','latex');
        title(['Differentiator reconstruction',newline,'automatic']);
        ax = gca;
        ax.FontSize = 25;
    end




    %% Plot attractor reconstruction error
    custom_map = [
            0.0 0.0 0.8;   % dark blue (low error)
            0.0 0.5 1.0;   % medium blue
            0.0 0.8 1.0;   % light blue
            1.0 0.6 0.0;   % orange
            1.0 0.0 0.0;   % red
            0.6 0.0 0.0    % dark red (high error)
        ];

    % Component-wise normalization
    nOrders = min(size(diff_output,2),3);
    gt_cell = cell(1,nOrders);
    gt_cell{1} = zero_der;
    if nOrders >= 2
        gt_cell{2} = first_der;
    end
    if nOrders >= 3
        gt_cell{3} = second_der;
    end

    sigma = zeros(1,nOrders);
    for k = 1:nOrders
        sigma(k) = std(gt_cell{k});
    end
    sigma(sigma == 0) = 1; % in case a ground-truth coordinate is constant throughout the signal, avoids division by zero

    % 2D embedding
    if size(diff_output,2) >= 2
        figID = figID + 1;
        figure(figID);
        % Ground-truth embedding
        y_gt_2d = [zero_der, first_der];
        % Errors for both methods
        y_guess_2d = [diff_output_guess(:,1), diff_output_guess(:,2)];
        y_auto_2d  = [diff_output(:,1),       diff_output(:,2)];

        % Component-wise normalized errors
        sigma_2d = sigma(1:2);
        err_guess_2d = sqrt(sum(((y_guess_2d - y_gt_2d)./sigma_2d).^2, 2));
        err_auto_2d  = sqrt(sum(((y_auto_2d  - y_gt_2d)./sigma_2d).^2, 2));

        % Common color scale
        cmax_2d = max([err_guess_2d; err_auto_2d]);
        if cmax_2d == 0
            cmax_2d = 1;
        end

        % Guess method
        subplot(1,2,1);
        ax_guess_2d = gca;
        scatter(ax_guess_2d, diff_output_guess(:,1), diff_output_guess(:,2), ...
            20, err_guess_2d, 'filled'); 
        xlabel('$x_{1,D}$','interpreter','latex');
        ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        title(['Guess parameters',newline,'reconstruction error']);
        ax_guess_2d.FontSize = 25;
        colormap(gcf, custom_map);
        clim(ax_guess_2d, [0 cmax_2d]);

        % Automatic method
        subplot(1,2,2);
        ax_auto_2d = gca;
        scatter(ax_auto_2d, diff_output(:,1), diff_output(:,2), ...
            20, err_auto_2d, 'filled'); 
        xlabel('$x_{1,D}$','interpreter','latex');
        ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        title(['Automatic',newline,'reconstruction error']);
        ax_auto_2d.FontSize = 25;
        colormap(gcf, custom_map);
        clim(ax_auto_2d, [0 cmax_2d]);

        % Common colorbar
        cb = colorbar(ax_auto_2d);
        cb.Label.String = '';
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize = 20;

        % Error histograms
        figID = figID + 1;
        figure(figID);
        ax_hist_2d = axes;
        hold(ax_hist_2d, 'on');

        [~, hist_edges_2d] = histcounts( ...
            [err_guess_2d; err_auto_2d], 'BinMethod', 'auto');

        histogram(ax_hist_2d, err_guess_2d, hist_edges_2d, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.00 0.35 0.85], ...
            'FaceAlpha', 0.45, ...
            'EdgeColor', 'none');

        histogram(ax_hist_2d, err_auto_2d, hist_edges_2d, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.75 0.00 0.00], ...
            'FaceAlpha', 0.45, ...
            'EdgeColor', 'none');

        xlabel(ax_hist_2d, 'Component-wise normalized L2 error', ...
            'Interpreter', 'latex');
        ylabel(ax_hist_2d, 'Probability', ...
            'Interpreter', 'latex');
        legend(ax_hist_2d, {'Guess parameters','Automatic'}, ...
            'Interpreter', 'latex', 'Location', 'best');

        ax_hist_2d.FontSize = 20;
        grid(ax_hist_2d, 'on');
        box(ax_hist_2d, 'on');

        stats_guess_2d = sprintf( ...
            'Guess\nMean: %.3f\nMedian: %.3f\nVariance: %.3f', ...
            mean(err_guess_2d), median(err_guess_2d), var(err_guess_2d));

        stats_auto_2d = sprintf( ...
            'Automatic\nMean: %.3f\nMedian: %.3f\nVariance: %.3f', ...
            mean(err_auto_2d), median(err_auto_2d), var(err_auto_2d));

        text(ax_hist_2d, 0.68, 0.95, stats_guess_2d, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [0.85 0.92 1.00], ...
            'EdgeColor', [0.00 0.35 0.85], ...
            'LineWidth', 1.2, ...
            'FontSize', 15, ...
            'Interpreter', 'none');

        text(ax_hist_2d, 0.68, 0.62, stats_auto_2d, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [1.00 0.88 0.88], ...
            'EdgeColor', [0.75 0.00 0.00], ...
            'LineWidth', 1.2, ...
            'FontSize', 15, ...
            'Interpreter', 'none');

            title('2D embedding');
    end

    % 3D embedding
    if size(diff_output,2) >= 3
        figID = figID + 1;
        figure(figID);
        % Ground-truth embedding
        y_gt_3d = [zero_der, first_der, second_der];
        % Errors for both methods
        y_guess_3d = [diff_output_guess(:,1), diff_output_guess(:,2), diff_output_guess(:,3)];
        y_auto_3d  = [diff_output(:,1),       diff_output(:,2),       diff_output(:,3)];

        % Component-wise normalized errors
        sigma_3d = sigma(1:3);
        err_guess_3d = sqrt(sum(((y_guess_3d - y_gt_3d)./sigma_3d).^2, 2));
        err_auto_3d  = sqrt(sum(((y_auto_3d  - y_gt_3d)./sigma_3d).^2, 2));

        % Normalized errors
        c_guess_3d = err_guess_3d;
        c_auto_3d  = err_auto_3d;

        % Common color scale
        cmax_3d = max([c_guess_3d; c_auto_3d]);
        if cmax_3d == 0
            cmax_3d = 1;
        end

        % Guess method (left)
        subplot(1,2,1);
        ax_guess_3d = gca;
        scatter3(ax_guess_3d, diff_output_guess(:,1), diff_output_guess(:,2), ...
            diff_output_guess(:,3), 20, c_guess_3d, 'filled');  % increased marker size
        xlabel('$x_{1,D}$','interpreter','latex');
        ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        zlabel('$\ddot{x}_{1,D}$','interpreter','latex');
        title(['Guess parameters',newline,'reconstruction error']);
        ax_guess_3d.FontSize = 25;
        colormap(gcf, custom_map);
        clim(ax_guess_3d, [0 cmax_3d]);

        % Automatic method (right)
        subplot(1,2,2);
        ax_auto_3d = gca;
        scatter3(ax_auto_3d, diff_output(:,1), diff_output(:,2), ...
            diff_output(:,3), 20, c_auto_3d, 'filled');  % increased marker size
        xlabel('$x_{1,D}$','interpreter','latex');
        ylabel('$\dot{x}_{1,D}$','interpreter','latex');
        zlabel('$\ddot{x}_{1,D}$','interpreter','latex');
        title(['Automatic',newline,'reconstruction error']);
        ax_auto_3d.FontSize = 25;
        colormap(gcf, custom_map);
        clim(ax_auto_3d, [0 cmax_3d]);

        % Common colorbar
        cb = colorbar(ax_auto_3d);
        cb.Label.String = '';
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize = 20;

        % Error histograms
        figID = figID + 1;
        figure(figID);
        ax_hist_3d = axes;
        hold(ax_hist_3d, 'on');

        [~, hist_edges_3d] = histcounts( ...
            [c_guess_3d; c_auto_3d], 'BinMethod', 'auto');

        histogram(ax_hist_3d, c_guess_3d, hist_edges_3d, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.00 0.35 0.85], ...
            'FaceAlpha', 0.45, ...
            'EdgeColor', 'none');

        histogram(ax_hist_3d, c_auto_3d, hist_edges_3d, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.75 0.00 0.00], ...
            'FaceAlpha', 0.45, ...
            'EdgeColor', 'none');

        xlabel(ax_hist_3d, 'Component-wise normalized L2 error', ...
            'Interpreter', 'latex');
        ylabel(ax_hist_3d, 'Probability', ...
            'Interpreter', 'latex');
        legend(ax_hist_3d, {'Guess parameters','Automatic'}, ...
            'Interpreter', 'latex', 'Location', 'best');

        ax_hist_3d.FontSize = 20;
        grid(ax_hist_3d, 'on');
        box(ax_hist_3d, 'on');

        stats_guess_3d = sprintf( ...
            'Guess\nMean: %.3f\nMedian: %.3f\nVariance: %.3f', ...
            mean(c_guess_3d), median(c_guess_3d), var(c_guess_3d));

        stats_auto_3d = sprintf( ...
            'Automatic\nMean: %.3f\nMedian: %.3f\nVariance: %.3f', ...
            mean(c_auto_3d), median(c_auto_3d), var(c_auto_3d));

        text(ax_hist_3d, 0.68, 0.95, stats_guess_3d, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [0.85 0.92 1.00], ...
            'EdgeColor', [0.00 0.35 0.85], ...
            'LineWidth', 1.2, ...
            'FontSize', 15, ...
            'Interpreter', 'none');

        text(ax_hist_3d, 0.68, 0.62, stats_auto_3d, ...
            'Units', 'normalized', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', [1.00 0.88 0.88], ...
            'EdgeColor', [0.75 0.00 0.00], ...
            'LineWidth', 1.2, ...
            'FontSize', 15, ...
            'Interpreter', 'none');

        title('3D embedding');
    end
    
end