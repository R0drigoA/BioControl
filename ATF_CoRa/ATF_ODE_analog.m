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


%figure; 
%plot(t, y(:, 2), '-')
%legend({"Y"})

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
    mY = 0.125;
    gV = 1;
    mV = 0.125;
    
     %Perturbation
    if (t>1000)
       mY = 10;       
    end



    %if (t>800)
    %    mY = 0.125;
    %end

    % Species:
    W = y(1);
    Y = y(2);
    V = y(3);
    U = y(4);
    C = y(5);

    % ODEs:   
    dWdt = mW - (g+gW)*W - np*U*W + (n0+gU)*C;
    dYdt = mY*W - (g+gY)*Y;
    dVdt = mY*((0.1+0.0008*C)/(0.0008+0.0375*U)) - (g+gV)*V;
    dUdt = mU*V - (g+gU)*U - np*U*W + (n0+gW)*C;
    dCdt = np*U*W - (g+n0+nm+gU+gW)*C;

    dydt = [dWdt; dYdt; dVdt; dUdt; dCdt];
end


