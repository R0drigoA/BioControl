%% ATF 1 CoRa 
% Code to calculate CoRa value for an ATF motif with mY = 0.125

clear; clc;

%% Simulate ODE dynamics to find the SS of the feedback system:
y0 = [0, 0, 0, 0]; % Initial conditions 
tspan = linspace(0, 50000, 5000);    
[t, y] = ode45(@f1, tspan, y0);

% Plot ODE results
figure; 
plot(t(500:5000), round(y(500:5000, 2),4), '-')
legend({"Y"})

SS = y(2499,:);     % array of steady states before perturbation
SS_p = y(5000,:);   % array of steady states after perturbation 

global W_ss;
W_ss = SS(1);       % SS of W before perturbation

Yss = SS(2);        % SS of Y before the perturbation
Yss_p = SS_p(2);    % SS of Y after the perturbation



%% Simulate ODE dynamics to find the SS of the analog system without feedback:

yp = [SS(1), SS(2), SS(2), SS(3), SS(4)]; % Initial conditions 
tspan = linspace(0, 50000, 5000);
[t, y] = ode45(@f2, tspan, yp);

%Plot ODE results.
figure; 
plot(t, round(y(:, 2),2), '-', t, round(y(:, 3),2), '-')
legend({"Y", "V"})

Yss_nf = Yss;           % SS before perturbation (is the same as Y_ss)
Yss_p_nf = y(5000,2);   % SS after perturation


%% CoRa value
CoRa = log10(Yss_p/Yss) / log10(Yss_p_nf/Yss_nf);

%% ODE for SS
function dydt = f1(t, y)
    % Kinetic parameters:  
    g = 0.0004;    % global dilution
    gU = 0.0004;   % U degradation rate
    gW = 0.0004;   % W degradation rate
    mU = 0.125;    % U synthesis rate
    mW = 0.1;      % W synthesis rate
    n0 = 0.0004;   % the spontaneous unbinding rate of U and W
    np = 0.0375;   % binding rate of U and W
    nm = 0.5;      % codegradation rate of U, W , in the complex form C
    
    gY = 1;        % Y degradation rate
    mY = 0.125;    % Y synthesis rate

    %Perturbation
    if (t>25000)
        % mY = 0.125*1.05 = 0.13125 
        mY = 0.13125;
    end

    % Species:
    W = y(1);
    Y = y(2);
    U = y(3);
    C = y(4);
    % ODEs:    
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dUdt = mU*Y - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dUdt; dCdt];
end


function dydt = f2(t, y)
    % Kinetic parameters:
    g = 0.0004;    % global dilution
    gU = 0.0004;   % U degradation rate
    gW = 0.0004;   % W degradation rate
    mU = 0.125;    % U synthesis rate
    mW = 0.1;      % W synthesis rate
    n0 = 0.0004;   % the spontaneous unbinding rate of U and W
    np = 0.0375;   % binding rate of U and W
    nm = 0.5;      % codegradation rate of U, W , in the complex form C
    
    gY = 1;        % Y degradation rate
    gV = 1;        % V degradation rate (same as gY)
    mY = 0.125;    % Y synthesis rate
    mV = 0.125;    % V synthesis rate (same as mY)
   
    % Perturbation
    if (t>25000)
        % mY = 0.125*1.05 = 0.13125 
        mY = 0.13125;
    end

    % Species:
    W = y(1);
    Y = y(2);
    V = y(3);
    U = y(4);
    C = y(5);

    % W_ss from the previous simulation
    global W_ss;

    % ODEs:   
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dVdt = mV*(W_ss) - (g+gV)*V;
    dUdt = mU*V - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dVdt; dUdt; dCdt];
end
