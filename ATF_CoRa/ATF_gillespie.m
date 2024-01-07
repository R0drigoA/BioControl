clear; clc;

% Condiciones iniciales
%    [W, Y, U, C] 
y0 = [0, 0, 0, 0]; 
total_time = 2000;

[t, W, Y, U, C] = gillespie_algorithm(y0, total_time);

% Plot the results for each species separately
figure;
plot(t, [W; Y; U; C], 'LineWidth', 0.25);
xlabel("Tiempo")
ylabel("Concentracion")
legend({"W", "Y", "U", "C"})

%figure;
%plot(t, [Y], 'LineWidth', 0.25);
%xlabel("Tiempo")
%ylabel("Concentracion")
%legend({"Y"})

function [t, W, Y, U, C] = gillespie_algorithm(y0, total_time)
    % Kinetic parameters
    mW = 0.1;
    g = 0.0004;
    gU = 0.0004;
    gW = 0.0004;
    mU = 0.125;
    n0 = 0.0004;
    np = 0.0375;
    nm = 0.5;
    gY = 1;
    mY = 0.125;


% Initialize variables
    t = 0;
    W = y0(1);
    Y = y0(2);
    U = y0(3);
    C = y0(4);

    time_points = [t];
    W_history = W;
    Y_history = Y;
    U_history = U;
    C_history = C;

    while t < total_time

        %Perturbation
        if (t>1000)
            mY = 0.5;
        end

        % Calculate reaction rates
        rs1 = mW + (n0 + gU) * C;
        rd1 = (g+ gW)*W + np * U * W;
        rs2 = mY * W;
        rd2 = (g +gY) * Y;
        rs3 = mU * Y + (n0 + gW) * C;
        rd3 = (g+ gU)*U + np * U * W;
        rs4 = np * U * W;
        rd4 = (g + n0 + nm + gU + gW) * C;
        

        % Calculate total reaction rate
        total_rate = rs1 + rd1 + rs2 + rd2 + rs3 + rd3 + rs4 + rd4;

        if total_rate == 0
            break; % No more reactions can occur
        end

        % Generate two random numbers
        r1 = rand();
        r2 = rand();

        % Calculate time until next reaction
        tau = (1 / total_rate) * log(1 / r1);

        % Determine which reaction occurs
        if r2 <= rs1 / total_rate
            W = W + 1;
        elseif r2 <= (rs1 + rd1) / total_rate
            W = W - 1;
        elseif r2 <= (rs1 + rd1 + rs2) / total_rate
            Y = Y + 1;
        elseif r2 <= (rs1 + rd1 + rs2 + rd2) / total_rate
            Y = Y - 1;
        elseif r2 <= (rs1 + rd1 + rs2 + rd2 + rs3) / total_rate
            U = U + 1;
        elseif r2 <= (rs1 + rd1 + rs2 + rd2 + rs3 + rd3) / total_rate
            U = U - 1;
        elseif r2 <= (rs1 + rd1 + rs2 + rd2 + rs3 + rd3 + rs4) / total_rate
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
        U_history = [U_history, U];
        C_history = [C_history, C];
    end

    % Return the time points and species counts
    t = time_points;
    W = W_history;
    Y = Y_history;
    U = U_history;
    C = C_history;
end


