%% Prueba para git
clear; clc;

% Simulate ODE dynamics:
y0 = [0, 0, 0]; % Initial conditions 
tspan = linspace(0, 500, 100);
[t, y] = ode45(@f, tspan, y0);

% Plot ODE results
figure;
plot(t, y(:, 1), '-', t, y(:, 2), '-', t, y(:, 3), '-')
xlabel("Tiempo")
ylabel("Concentracion")
legend({"A", "B", "C"})

function dydt = f(t, y)
    % Kinetic parameters:
    gA = 0.0004;
    gB = 0.0004;
    mA = 0.125;
    mB = 0.1;
    n0 = 0.0004;
    np = 0.0375;
    nm = 0.5;

    % Perturbation
    if (t > 200)
        mA = 0.5;
    end

    % Species:
    A = y(1);
    B = y(2);
    C = y(3);

    % ODEs:
    dadt = mA - gA * A - np * A * B + n0 * C;
    dbdt = mB * A - gB * B - np * A * B + n0 * C;
    dcdt = np * A * B - n0 * C - nm * C;

    dydt = [dadt; dbdt; dcdt];
end

