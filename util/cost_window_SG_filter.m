function cost_value = cost_window_SG_filter(window_length,m,lambda,t,chattering,diff_output_opt,norm_factor_chat,norm_factor_dist)
% This function computes the SG cost for an SG window size window_length.

        % Convert window into number of data points
        [~,window_length] = min(abs((t-t(1))-window_length));
        if mod(window_length,2) == 0
            window_length = window_length + 1;
        end
        % Apply SG filter to the chattering component and optimal Differentiator output
        smoothed_chat = sgolayfilt(chattering,m,window_length);
        smoothed_diff = sgolayfilt(diff_output_opt,m,window_length);

        % Compute cost components
        cost_chat = 1/(length(smoothed_chat)-1)*sum(abs(diff(smoothed_chat)));
        cost_dist = 1/var(smoothed_diff);
   
        % Trade-off cost window lenght computation
        cost_chat = 1/norm_factor_chat*cost_chat;
        cost_dist = 1/norm_factor_dist*cost_dist;
    
        cost_value = (1-lambda)*cost_chat + lambda*cost_dist;

    end