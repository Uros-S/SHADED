function L_opt = autotune_Differentiator(signal,t,diff_parameters,complexity_settings,iteration)
% This function automatically computes the Differentiator gain L based on
% the Differentiator cost.
    
    t_init = tic;

    % Differentiator parameters
    n_d = diff_parameters{2}(iteration);
    n_f = diff_parameters{3}(iteration);
    q = diff_parameters{4}(iteration);

    if isnan(n_f)
        n_f = 0;
    end

    % Compute guess L_0 with finite difference derivatives
    tmp = signal;
    for i=1:n_d+1
        tmp = finite_diff_derivative(tmp,t,1);
    end
    L_0_guess = max(abs(tmp));

    L_0_max = 10*L_0_guess;

    cost = @(x) cost_differentiator(signal,t,n_d,n_f,x,q,L_0_max);

    % Find optimal L_0 value that minimises the residual variance with simulated annealing
    max_iterations_diff = complexity_settings{1}(iteration);
    if isnan(max_iterations_diff)
        max_iterations_diff = 150;
    end
    tolerance_diff = complexity_settings{2}(iteration);
    if isnan(tolerance_diff)
        tolerance_diff = 1e-4;
    end
    SA_options = optimoptions('simulannealbnd','MaxIterations',max_iterations_diff,'FunctionTolerance',tolerance_diff,'Display','off');
    if complexity_settings{7}(iteration)
        rng(123,'twister')   % for repeatability
    end
    x_opt = simulannealbnd(cost,log10(L_0_guess)/log10(L_0_max),0,1,SA_options);
    L_opt = L_0_max^x_opt;

    t_fin = toc(t_init);
    disp(['Differentiator tuning: ',num2str(round(t_fin,2,'significant')),' [s]']);
    disp(['Estimated L value is ',num2str(L_opt,'%.2g')]);
    
end