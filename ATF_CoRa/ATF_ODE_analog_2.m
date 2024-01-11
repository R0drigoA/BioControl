%% ATF 1 analog
clear; clc;

%% Simulate ODE dynamics to find the SS of the feedback system:
y0 = [0, 0, 0, 0]; % Initial conditions 
tspan = linspace(0, 50000, 5000);    
[t, y] = ode45(@f1, tspan, y0);

% Plot ODE results
figure; 
plot(t, y(:, 1), '-', t, y(:, 2), '-',  t, y(:, 3), '-', t, y(:, 4), '-')
legend({"W", "Y", "U", "C"})

SS = y(2499,:);
SS_p = y(5000,:);

Yss = SS(2);      % steady state value before the perturbation
Yss_p = SS_p(2);    % steady state value after the perturbation

%% Simulate ODE dynamics to find the SS of the analog system without feedback:
yp = [SS(1), SS(2), SS(2), SS(3), SS(4)]; % Initial conditions 
tspan = linspace(0, 50000, 5000);
[t, y] = ode45(@f2, tspan, yp);

%Plot ODE results.
figure; 
plot(t, y(:, 1), '-', t, y(:, 2), '-', t, y(:, 3), '-', t, y(:, 4), '-', t, y(:, 5), '-')
legend({"W", "Y", "V", "U", "C"})

figure; 
plot(t, y(:, 2), '-', t, y(:, 3), '-')
legend({"Y", "V"})

Yss_nf = Yss;
Yss_p_nf = y(5000,2);

CoRa = log10(Yss_p/Yss) / log10(Yss_p_nf/Yss_nf);

%% ODE for SS
function dydt = f1(t, y)
    % Kinetic parameters:
    g = 0.0004;
    gU = 0.0004;
    gW = 0.0004;
    mU = 0.125;
    mW = 0.1;
    n0 = 0.0004;
    np = 0.0375;
    nm = 0.5;
    gY = 1;   
    mY = 0.00125;

    %Perturbation
    if (t>25000)
       mY = 0.005;
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
    g = 0.0004;
    gU = 0.0004;
    gW = 0.0004;
    mU = 0.125;
    mW = 0.1;
    n0 = 0.0004;
    np = 0.0375;
    nm = 0.5;
    gY = 1;
    mY = 0.5;
    mV = 0.5;
    
    % Species:
    W = y(1);
    Y = y(2);
    V = y(3);
    U = y(4);
    C = y(5);
    W_ss = 104.5854882610933;

    % ODEs:   
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dVdt = mV*(W_ss) - (g+gY)*V;
    dUdt = mU*V - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dVdt; dUdt; dCdt];
end
