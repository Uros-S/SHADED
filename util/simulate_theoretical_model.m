function [t,s,ground_truth,model_param] = simulate_theoretical_model(model)
% This function simulates the model specified by the integer variable "model"

    t_init = 0;

    switch model

        case 1      % Lorenz
            sigma = 10; rho = 28; beta = 8/3;
            model_param = [sigma,rho,beta];
            dt = 1e-4;
            t_fin = 1020;
            t_trans = 1000;  
            x_0 = [1;1;1];
            options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

        case 2      % Hindmarsh-Rose
            a = 1; c = 1; d = 5; r = 0.01; sdyn = 4; x_R = -(1+sqrt(5))/2; I = 3; b = 2.7;
            model_param = [a,b,c,d,r,sdyn,I,x_R];
            dt = 1e-3;
            t_fin = 1900;
            t_trans = 1500;  
            x_0 = [0.65; 0.55; 0.45];
            options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

        case 3      % Fitzhugh-Nagumo
            a = 0.7; b = 0.8; I_ext = 0.5; tau = 12.5;
            model_param = [a,b,I_ext,tau];
            dt = 1e-3;
            t_fin = 190; 
            t_trans = 40;
            x_0 = [0,0];
            options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,2));

        case 4      % Jansen-Rit
            C = 135; A = 3.25; B = 22; a = 100; b = 50; V_0 = 6; nu_max = 5; r = 0.56; p = 200;
            model_param = [C,A,B,a,b,V_0,nu_max,r,p];
            dt = 1e-5; %1e-5;
            t_fin = 3; %3;
            t_trans = 1.75;
            x_0 = [0,0,0,0,0,0]; 
            options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,6));

        otherwise
            error('Model to be simulated not specified correctly!');
    end


    t_span = t_init:dt:t_fin;

    switch model
        case 1
            [t,s] = ode45(@(t,s)lorenz(t,s,sigma,rho,beta),t_span,x_0,options);
        case 2
            [t,s] = ode45(@(t,s)neuronHR(t,s,a,b,c,d,I,r,sdyn,x_R),t_span,x_0,options);
        case 3
            [t,s] = ode45(@(t,s)Fitzhugh_Nagumo(t,s,a,b,I_ext,tau),t_span,x_0,options);
        case 4
            [t,s] = ode45(@(t,s)Jansen_Rit(t,s,model_param),t_span,x_0,options);
    end
        
    cut_off_indx = find(t > t_trans,1);
    t= t(cut_off_indx:end);
    s = s(cut_off_indx:end,:);

    % Signal of interest
    switch model
        case {1,3}
            ground_truth = s(:,1);
        case 2
            ground_truth = s(:,1);
        case 4
            ground_truth = s(:,2)-s(:,3);
        case 5
            ground_truth = s(:,1);
    end

end