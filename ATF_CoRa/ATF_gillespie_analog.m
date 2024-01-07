clear; clc;

% Condiciones iniciales
%    [W, Y, V, U, C] 
y0 = [0, 0, 0, 0, 0]; 
total_time = 2000;

[t, W, Y, V, U, C] = gillespie_algorithm(y0, total_time);

% Plot the results for each species separately
figure;
plot(t, [W; Y; V; U; C], 'LineWidth', 0.25);
xlabel("Tiempo")
ylabel("Concentracion")
legend({"W", "Y", "V", "U", "C"})

% figure;
% plot(t, [Y], 'LineWidth', 0.25);
% xlabel("Tiempo")
% ylabel("Concentracion")
% legend({"Y"})

function [t, W, Y, V, U, C] = gillespie_algorithm(y0, total_time)
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


% Initialize variables
    t = 0;
    W = y0(1);
    Y = y0(2);
    V = y0(3);
    U = y0(4);
    C = y0(5);

    time_points = t;
    W_history = W;
    Y_history = Y;
    V_history = V;
    U_history = U;
    C_history = C;

    while t < total_time

        %Perturbation
        if (t>1000)
            mY = 0.5;
        end

        % Calculate reaction rates
        rates(1) = mW + (n0+gU)*C;    % rs W
        rates(2) = (g+gW)*W + np*U*W; % rd W
        rates(3) = mY*W;              % rs Y
        rates(4) = (g+gY)*Y;          % rd Y
        rates(5) = mV*W;              % rs V 
        rates(6) = (g+gV)*V;          % rd V
        rates(7) = mU*V + (n0+gW)*C;  % rs U
        rates(8) = (g+gU)*U + np*U*W; % rd U
        rates(9) = np*U*W;            % rs C
        rates(10) = (g+n0+nm+gU+gW)*C;% rd C
        

        % Calculate total reaction rate
        total_rate = sum(rates);

        if total_rate == 0
            break; % No more reactions can occur
        end

        % Generate two random numbers
        r1 = rand();
        r2 = rand();

        % Calculate time until next reaction
        tau = (1 / total_rate) * log(1 / r1);

        % Determine which reaction occurs
        if r2 <= rates(1) / total_rate
            W = W + 1;
        elseif r2 <= sum(rates(1:2)) / total_rate
            W = W - 1;
        elseif r2 <= sum(rates(1:3)) / total_rate
            Y = Y + 1;
        elseif r2 <= sum(rates(1:4)) / total_rate
            Y = Y - 1;
        elseif r2 <= sum(rates(1:5)) / total_rate
            V = V + 1;
        elseif r2 <= sum(rates(1:6)) / total_rate
            V = V - 1;
        elseif r2 <= sum(rates(1:7)) / total_rate
            U = U + 1;
        elseif r2 <= sum(rates(1:8)) / total_rate
            U = U - 1;
        elseif r2 <= sum(rates(1:9)) / total_rate
            C = C + 1; 
        else 
            C = C - 1; 
        end

        % Update time and species counts
        t = t + tau;

        % Store time and species values
        time_points = [time_points, t];
        W_history = [W_history, W];
        Y_history = [Y_history, Y];
        V_history = [V_history, V];
        U_history = [U_history, U];
        C_history = [C_history, C];
    end

    % Return the time points and species counts
    t = time_points;
    W = W_history;
    Y = Y_history;
    V = V_history;
    U = U_history;
    C = C_history;
end


