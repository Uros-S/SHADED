% Uros Sutulovic, 09/2026

clear; close all; clc;
addpath([pwd,'/systems/']) 
addpath([pwd,'/util/'])
figID = 0;

%% Selection of models, noise and result to be obtained

model = 1;                  % 1 = Lorenz
                            % 2 = Hindmarsh-Rose
                            % 3 = Fitzhugh-Nagumo
                            % 4 = Jansen-Rit

modality = 1;               % 1 = estimate derivatives with selected parameters
                            % 2 = investigate parameter L
                            % 3 = investigate SG filter window size for a given L
                            % 4 = sensitivity analysis of methodology parameters

noise_type = 1;             % 0 = no noise
                            % 1 = white Gaussian

additive_noise = 1;         % 0 = multiplicative noise
                            % 1 = additive noise

noise_var = 1;              % Variance of white Gaussian noise if noise_type == 1

n_derivatives = 2;          % Derivative order to estimate (only for modality 1)

persistence = true;         % Enable or disable frequency-persistence criterion

%% Parameters setting for methodology with guess parameter values
% No staircase architecture, the same parameters are applied for all the
% derivative orders

% Differentiator guess parameters
L_guess = 300000; 
n_d_guess = 2;
n_f_guess = 0;                  
q_guess = 0;                 % NaN to disable discrete L-adaptation

% SG filter guess parameters
m_guess = 2;
fl_guess = 0;   % in number of samples 

%% Parameter settings for Differentiator, Savitzky-Golay filter and computational complexity
% Select NaN for the parameter to be estimated/set to default settings. A
% scalar parameter selection will apply that value to all the levels of the
% staircase up to the derivative order to estimate n_derivatives, for
% customized choices select a vector with entries corresponding to the
% derivative order starting from the zeroth-order derivative (if the vector 
% has length less than n_derivatives, the remaining entries will be filled
% with NaN)

% Differentiator parameters
n_d = NaN;
L = NaN;  

n_f = NaN;                  
q = 0;                 % NaN to disable discrete L-adaptation
trans_diff = NaN;      % Transitory time to be discarded (in time units)

diff_parameters = {L,n_d,n_f,q,trans_diff};

% Savitzky-Golay filtering of Differentiator output
m = NaN;
fl = NaN;
SG_n_points_min = NaN;   % in number of samples   
SG_n_points_max = NaN;   % in number of samples 
SG_iterations = NaN;     
lambda = NaN;
n_windows = NaN;
overlap = NaN;
epsilon = NaN;
SG_parameters = {m,fl,SG_n_points_min,SG_n_points_max,SG_iterations,lambda,n_windows,overlap,epsilon};

% Computational complexity settings
max_iterations_diff = NaN; 
tolerance_diff = NaN;
max_iterations_SG = NaN;
tolerance_SG = NaN;
SG_cost = true;
use_seed = true;
complexity_settings = {max_iterations_diff,tolerance_diff,max_iterations_SG,tolerance_SG,persistence,SG_cost,use_seed};

%% Simulate theoretical model and add noise
[t_original,s,ground_truth,model_param] = simulate_theoretical_model(model);

% Data corruption
measurement = synthetic_data_corruption(ground_truth,noise_type,noise_var,additive_noise);

switch modality
    %% Differentiator plus Savitzky-Golay filtering attractor reconstruction
    case 1
        % Derivative estimation by guessing parameters
        if n_d_guess < n_derivatives
            n_d_guess = n_derivatives;
        end

        t_start_D_guess = tic;
        diff_signals_guess = differentiator(measurement,t_original,n_d_guess,n_f_guess,L_guess,q_guess,0);
        if fl_guess~=0
            diff_signals_guess = sgolayfilt(diff_signals_guess,m_guess,fl_guess);
        end

        t_end_D_guess = toc(t_start_D_guess);
        disp(['Computation time Differentiator plus Savitzky-Golay filter with guessed parameters: ',num2str(round(t_end_D_guess,2,'significant')),' [s]']);
        disp(newline)

        % Estimation of derivatives with automatic method
        [t,diff_signals] = automatic_derivative(t_original,measurement,n_derivatives,diff_parameters,SG_parameters,complexity_settings);

        if length(t) < length(t_original)
            indx1 = find(t_original==t(1),1);
            indx2 = find(t_original==t(end),1);
            s = s(indx1:indx2,:);
            ground_truth = ground_truth(indx1:indx2,:);
            measurement = measurement(indx1:indx2,:);
            diff_signals_guess = diff_signals_guess(indx1:indx2,:);
        end
        
        % Plot results
        figID = plot_signals_theoretical_model(model,t,s,measurement,diff_signals,diff_signals_guess,model_param,SG_parameters,figID);

    %% Analysis of the Differentiator gain inference
    case 2

        exp_min = 0;
        if model == 1      % 1 = Lorenz
            exp_max = 20;
        elseif model == 2  % 2 = Hindmarsh-Rose
            exp_max = 20;
        elseif model == 3  % 3 = Fitzhugh-Nagumo
            exp_max = 12;
        elseif model == 4  % 4 = Jansen-Rit
            exp_max = 25;
        else
            error('Model selected not valid!');
        end
        N_L = 1000;
        L_settings = [exp_min,exp_max,N_L];

        if isnan(n_d)
            error('Specify Differentiator order!');
        end

        [results_L,figID] = sweep_L(t_original,ground_truth,measurement,diff_parameters,L_settings,figID);

    %% Analysis of SG filter window size inference
    case 3
        % Effect of SG filter window size
        length_min = m+1;
        if model == 1
            length_max = 1;           % in time units
        elseif model == 2
            length_max = 2;           % in time units
        elseif model == 3
            length_max = 8;           % in time units
        elseif model == 4
            length_max = 0.125;       % in time units
        else
            error('Model selected not valid!');
        end

        length_max = round(length_max/(t_original(2)-t_original(1)));   % convert length in time units to number of samples
        N_SG = 1000;

        fl_settings = [length_min,length_max,N_SG];

        diff_parameters{2} = 0;     % consider denoising case
        if isnan(L)
            error('Specify Differentiator gain for sensitivity analysis!');
        end

        [results_SG,SG_cost,freq_max_psd,figID] = sweep_SG_window(t_original,ground_truth,measurement,fl_settings,diff_parameters,SG_parameters,figID);

    %% Sensitivity analysis of parameters
    case 4

        param_number = 1;   % 1 = lambda
                            % 2 = epsilon
                            % 3 = w_max

        N_param = 20;   % how many data points to consider

        if param_number == 1
            param_min = 0.1;
            param_max = 0.9;
        elseif param_number == 2
            param_min = 0.001;
            param_max = 0.2;
        elseif param_number == 3
            param_min = 10000;
            param_max = 20001;
        else
            error('Parameter selected for sensitivity analysis not valid!');
        end

        param_settings = [param_min,param_max,N_param];

        diff_parameters{2} = 0;     % consider denoising case
        if isnan(L)
            error('Specify Differentiator gain for sensitivity analysis!');
        end

        [results_param,figID] = sweep_parameter(t_original,measurement,diff_parameters,SG_parameters,complexity_settings,param_settings,param_number,figID);

    otherwise
        error('Modality inserted not valid!');
end

