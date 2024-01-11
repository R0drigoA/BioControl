% hola
clear; clc;

% Simulate ODE dynamics:
y0 = [0, 0, 0, 0, 0]; % Initial conditions 
tspan = linspace(0, 2000, 100);
[t, y] = ode45(@f, tspan, y0);

% Plot ODE results
figure; 
plot(t, y(:, 1), '-', t, y(:, 2), '-', t, y(:, 3), '-', t, y(:, 4), '-', t, y(:, 5), '-')
legend({"W", "Y", "V", "U", "C"})

function dydt = f(t, y)
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
    mY = 0.05;

    % Species:
    W = y(1);
    Y = y(2);
    V = y(3);
    U = y(4);
    C = y(5);

     %Perturbation
    if (t>1000)
       mY = 0.08;
    end

    % ODEs:   
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dVdt = mY*((mW+(g+gU)*C)/(g+gW+np*U)) - (g+gY)*V;
    dUdt = mU*V - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dVdt; dUdt; dCdt];
end


