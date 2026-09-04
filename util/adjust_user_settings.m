function [diff_parameters,SG_parameters,complexity_settings] = adjust_user_settings(n_derivatives,diff_parameters,SG_parameters,complexity_settings)
% This function accomodates manually selected parameters in the proposed
% methodology for every stage of the staircase architecture needed.

    %% Differentiator parameters
    if length(diff_parameters{1}) > n_derivatives + 1
        error('Number of Differentiator paremeters for L_0 is more than the derivative order to be estimated!');
    elseif length(diff_parameters{1}) == 1
        diff_parameters{1} = ones(1,n_derivatives + 1)*diff_parameters{1};
    elseif length(diff_parameters{1}) < n_derivatives
        diff_parameters{1} = [diff_parameters{1},nan(1,n_derivatives-length(diff_parameters{1}) + 1)];
    end

    if length(diff_parameters{2}) > n_derivatives + 1
        error('Number of Differentiator paremeters for n_d is more than the derivative order to be estimated!');
    elseif length(diff_parameters{2}) == 1
        diff_parameters{2} = ones(1,n_derivatives + 1)*diff_parameters{2};
    elseif length(diff_parameters{2}) < n_derivatives
        diff_parameters{2} = [diff_parameters{2},nan(1,n_derivatives-length(diff_parameters{2}) + 1)];
    end

    if length(diff_parameters{3}) > n_derivatives + 1
        error('Number of Differentiator paremeters for n_f is more than the derivative order to be estimated!');
    elseif length(diff_parameters{3}) == 1
        diff_parameters{3} = ones(1,n_derivatives + 1)*diff_parameters{3};
    elseif length(diff_parameters{3}) < n_derivatives
        diff_parameters{3} = [diff_parameters{3},nan(1,n_derivatives-length(diff_parameters{3}) + 1)];
    end

    if length(diff_parameters{4}) > n_derivatives + 1
        error('Number of Differentiator paremeters for q is more than the derivative order to be estimated!');
    elseif length(diff_parameters{4}) == 1
        diff_parameters{4} = ones(1,n_derivatives + 1)*diff_parameters{4};
    elseif length(diff_parameters{4}) < n_derivatives
        diff_parameters{4} = [diff_parameters{4},nan(1,n_derivatives-length(diff_parameters{4}) + 1)];
    end

    if length(diff_parameters{5}) > n_derivatives + 1
        error('Number of Differentiator paremeters for trans_diff is more than the derivative order to be estimated!');
    elseif length(diff_parameters{5}) == 1
        diff_parameters{5} = ones(1,n_derivatives + 1)*diff_parameters{5};
    elseif length(diff_parameters{5}) < n_derivatives
        diff_parameters{5} = [diff_parameters{5},nan(1,n_derivatives-length(diff_parameters{5}) + 1)];
    end

    %% Savitzky-Golay filter parameters
    if length(SG_parameters{1}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for m is more than the derivative order to be estimated!');
    elseif length(SG_parameters{1}) == 1
        SG_parameters{1} = ones(1,n_derivatives + 1)*SG_parameters{1};
    elseif length(SG_parameters{1}) < n_derivatives
        SG_parameters{1} = [SG_parameters{1},nan(1,n_derivatives-length(SG_parameters{1}) + 1)];
    end

    if length(SG_parameters{2}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for fl is more than the derivative order to be estimated!');
    elseif length(SG_parameters{2}) == 1
        SG_parameters{2} = ones(1,n_derivatives + 1)*SG_parameters{2};
    elseif length(SG_parameters{2}) < n_derivatives
        SG_parameters{2} = [SG_parameters{2},nan(1,n_derivatives-length(SG_parameters{2}) + 1)];
    end

    if length(SG_parameters{3}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for SG_n_points_min is more than the derivative order to be estimated!');
    elseif length(SG_parameters{3}) == 1
        SG_parameters{3} = ones(1,n_derivatives + 1)*SG_parameters{3};
    elseif length(SG_parameters{3}) < n_derivatives
        SG_parameters{3} = [SG_parameters{3},nan(1,n_derivatives-length(SG_parameters{3}) + 1)];
    end

    if length(SG_parameters{4}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for SG_n_points_max is more than the derivative order to be estimated!');
    elseif length(SG_parameters{4}) == 1
        SG_parameters{4} = ones(1,n_derivatives + 1)*SG_parameters{4};
    elseif length(SG_parameters{4}) < n_derivatives
        SG_parameters{4} = [SG_parameters{4},nan(1,n_derivatives-length(SG_parameters{4}) + 1)];
    end

    if length(SG_parameters{5}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for delta_points is more than the derivative order to be estimated!');
    elseif length(SG_parameters{5}) == 1
        SG_parameters{5} = ones(1,n_derivatives + 1)*SG_parameters{5};
    elseif length(SG_parameters{5}) < n_derivatives
        SG_parameters{5} = [SG_parameters{5},nan(1,n_derivatives-length(SG_parameters{5}) + 1)];
    end

    if length(SG_parameters{6}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for lambda is more than the derivative order to be estimated!');
    elseif length(SG_parameters{6}) == 1
        SG_parameters{6} = ones(1,n_derivatives + 1)*SG_parameters{6};
    elseif length(SG_parameters{6}) < n_derivatives
        SG_parameters{6} = [SG_parameters{6},nan(1,n_derivatives-length(SG_parameters{6}) + 1)];
    end

    if length(SG_parameters{7}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for n_windows is more than the derivative order to be estimated!');
    elseif length(SG_parameters{7}) == 1
        SG_parameters{7} = ones(1,n_derivatives + 1)*SG_parameters{7};
    elseif length(SG_parameters{7}) < n_derivatives
        SG_parameters{7} = [SG_parameters{7},nan(1,n_derivatives-length(SG_parameters{7}) + 1)];
    end

    if length(SG_parameters{8}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for overlap is more than the derivative order to be estimated!');
    elseif length(SG_parameters{8}) == 1
        SG_parameters{8} = ones(1,n_derivatives + 1)*SG_parameters{8};
    elseif length(SG_parameters{8}) < n_derivatives
        SG_parameters{8} = [SG_parameters{8},nan(1,n_derivatives-length(SG_parameters{8}) + 1)];
    end

    if length(SG_parameters{9}) > n_derivatives + 1
        error('Number of Savitzky-Golay filter paremeters for epsilon is more than the derivative order to be estimated!');
    elseif length(SG_parameters{9}) == 1
        SG_parameters{9} = ones(1,n_derivatives + 1)*SG_parameters{9};
    elseif length(SG_parameters{9}) < n_derivatives
        SG_parameters{9} = [SG_parameters{9},nan(1,n_derivatives-length(SG_parameters{9}) + 1)];
    end

    %% Computational complexity settings
    if length(complexity_settings{1}) > n_derivatives + 1
        error('Number of computational complexity paremeters for max_iterations_diff is more than the derivative order to be estimated!');
    elseif length(complexity_settings{1}) == 1
        complexity_settings{1} = ones(1,n_derivatives + 1)*complexity_settings{1};
    elseif length(complexity_settings{1}) < n_derivatives
        complexity_settings{1} = [complexity_settings{1},nan(1,n_derivatives-length(complexity_settings{1}) + 1)];
    end

    if length(complexity_settings{2}) > n_derivatives + 1
        error('Number of computational complexity paremeters for tolerance_diff is more than the derivative order to be estimated!');
    elseif length(complexity_settings{2}) == 1
        complexity_settings{2} = ones(1,n_derivatives + 1)*complexity_settings{2};
    elseif length(complexity_settings{2}) < n_derivatives
        complexity_settings{2} = [complexity_settings{2},nan(1,n_derivatives-length(complexity_settings{2}) + 1)];
    end

    if length(complexity_settings{3}) > n_derivatives + 1
        error('Number of computational complexity paremeters for max_iterations_SG is more than the derivative order to be estimated!');
    elseif length(complexity_settings{3}) == 1
        complexity_settings{3} = ones(1,n_derivatives + 1)*complexity_settings{3};
    elseif length(complexity_settings{3}) < n_derivatives
        complexity_settings{3} = [complexity_settings{3},nan(1,n_derivatives-length(complexity_settings{3}) + 1)];
    end

    if length(complexity_settings{4}) > n_derivatives + 1
        error('Number of computational complexity paremeters for tolerance_SG is more than the derivative order to be estimated!');
    elseif length(complexity_settings{4}) == 1
        complexity_settings{4} = ones(1,n_derivatives + 1)*complexity_settings{4};
    elseif length(complexity_settings{4}) < n_derivatives
        complexity_settings{4} = [complexity_settings{4},nan(1,n_derivatives-length(complexity_settings{4}) + 1)];
    end

    if length(complexity_settings{5}) > n_derivatives + 1
        error('Number of computational complexity paremeters for persistence is more than the derivative order to be estimated!');
    elseif length(complexity_settings{5}) == 1
        complexity_settings{5} = ones(1,n_derivatives + 1)*complexity_settings{5};
    elseif length(complexity_settings{5}) < n_derivatives
        complexity_settings{5} = [complexity_settings{5},nan(1,n_derivatives-length(complexity_settings{5}) + 1)];
    end

    if length(complexity_settings{6}) > n_derivatives + 1
        error('Number of computational complexity paremeters for SG_cost is more than the derivative order to be estimated!');
    elseif length(complexity_settings{6}) == 1
        complexity_settings{6} = ones(1,n_derivatives + 1)*complexity_settings{6};
    elseif length(complexity_settings{6}) < n_derivatives
        complexity_settings{6} = [complexity_settings{6},nan(1,n_derivatives-length(complexity_settings{6}) + 1)];
    end

    if length(complexity_settings{7}) > n_derivatives + 1
        error('Number of computational complexity paremeters for use_seed is more than the derivative order to be estimated!');
    elseif length(complexity_settings{7}) == 1
        complexity_settings{7} = ones(1,n_derivatives + 1)*complexity_settings{7};
    elseif length(complexity_settings{7}) < n_derivatives
        complexity_settings{7} = [complexity_settings{7},nan(1,n_derivatives-length(complexity_settings{7}) + 1)];
    end

end