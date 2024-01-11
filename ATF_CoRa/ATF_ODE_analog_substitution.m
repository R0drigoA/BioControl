%% ATF analog 1
clear; clc;

%% Simulate ODE dynamics to find the SS:
y0 = [0, 0, 0, 0]; % Initial conditions 
tspan = linspace(0, 2000, 100);
[t, y] = ode45(@f1, tspan, y0);

% Plot ODE results
figure; 
plot(t, y(:, 1), '-', t, y(:, 2), '-',  t, y(:, 3), '-', t, y(:, 4), '-')
legend({"W", "Y", "U", "C"})

%figure; 
%plot(t, y(:, 2), '-')
%legend({"Y"})
W_ss = y(100,1);

%% Simulate ODE dynamics to find the SS of the analog system:
y_ss = [y(100,1), y(100,2), y(100,3), y(100,4)]; % Initial conditions 
tspan = linspace(0, 2000, 100);
[t, y] = ode45(@f2, tspan, y_ss);

% Plot ODE results.
figure; 
plot(t, y(:, 1), '-', t, y(:, 2), '-', t, y(:, 3), '-', t, y(:, 4), '-')
legend({"W", "Y", "U", "C"})

%figure; 
%plot(t, y(:, 2), '-')
%legend({"Y"})% using a new constant

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
    mY = 0.125;

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

    % Species:
    W = y(1);
    Y = y(2);
    U = y(3);
    C = y(4);

    % ODEs:   
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dUdt = mU*(2) - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dUdt; dCdt];
end





