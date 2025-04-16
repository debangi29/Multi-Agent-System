% MATLAB code for simulation Example 2 from the thesis
% "Fixed Time Consensus Analysis of Multi Agent Systems via Impulsive Control"
% Case: Impulsive strength greater than one (with and without necessary conditions)

clear all;
close all;

% Define the adjacency matrix A and Laplacian matrix L for example 2
A = [0 1 1 0;
     1 0 1 1;
     1 1 0 0;
     0 1 0 0];

L = [2 -1 -1 0;
    -1 3 -1 -1;
    -1 -1 2 0;
     0 -1 0 1];

% Parameters for the control protocol
c1 = 2;
c2 = 2;
c3 = 2;
alpha = 3/5;
gamma = 11/7;

% Initial states of the agents for Example 2
x0 = [5; -2; 6; -1];

% Simulation parameters
T = 1.0; % Total simulation time
dt = 0.001; % Time step
t = 0:dt:T;
N = length(t);
n = length(x0); % Number of agents

% Function to calculate sgn^alpha(x)
sgn_alpha = @(x) sign(x).*abs(x).^alpha;

% Function to calculate sgn^gamma(x)
sgn_gamma = @(x) sign(x).*abs(x).^gamma;

% Eigen analysis for the necessary condition
[V, D] = eig(L);
eigenvalues = diag(D);
eigenvalues = sort(eigenvalues);
lambda_2 = eigenvalues(2); % Second smallest eigenvalue
lambda_N = eigenvalues(end); % Largest eigenvalue

% Arrays to store the states for both cases
X_with_condition = zeros(n, N);
X_without_condition = zeros(n, N);

% Initial conditions
X_with_condition(:, 1) = x0;
X_without_condition(:, 1) = x0;

% Impulsive parameters
impulsive_strength = 1.2; % Impulsive strength > 1
critical_value = 2*lambda_2/(lambda_N^2);
fprintf('Critical condition for stability: impulsive_strength > %.4f\n', critical_value);

% Critical theta_2 based on necessary condition
critical_theta_2 = log(impulsive_strength)/c3;
fprintf('Critical value for theta_2: %.4f\n', critical_theta_2);

% Case 1: With necessary conditions satisfied
theta_2_good = 0.1;  % Chosen to be greater than critical_theta_2
fprintf('Case 1 (with necessary conditions): theta_2 = %.2f > %.4f\n', theta_2_good, critical_theta_2);

% Case 2: Without necessary conditions satisfied
theta_2_bad = 0.05; % Chosen to be less than critical_theta_2
fprintf('Case 2 (without necessary conditions): theta_2 = %.3f < %.4f\n', theta_2_bad, critical_theta_2);

% Define the impulsive times for both cases
impulsive_times_good = 0:theta_2_good:T;
impulsive_times_bad = 0:theta_2_bad:T;

% Find indices for impulsive times
is_impulse_time_good = false(1, N);
for i = 1:length(impulsive_times_good)
    [~, idx] = min(abs(t - impulsive_times_good(i)));
    is_impulse_time_good(idx) = true;
end

is_impulse_time_bad = false(1, N);
for i = 1:length(impulsive_times_bad)
    [~, idx] = min(abs(t - impulsive_times_bad(i)));
    is_impulse_time_bad(idx) = true;
end

% Simulation with necessary conditions satisfied
for k = 1:N-1
    % Calculate control input for each agent
    for i = 1:n
        u_i = 0;
        for j = 1:n
            if A(i,j) > 0
                diff_ij = X_with_condition(j, k) - X_with_condition(i, k);
                u_i = u_i + A(i,j) * (c1 * sgn_alpha(diff_ij) + ...
                                     c2 * sgn_gamma(diff_ij) + ...
                                     c3 * sign(diff_ij));
            end
        end
        
        % Update state with continuous dynamics
        X_with_condition(i, k+1) = X_with_condition(i, k) + dt * u_i;
    end
    
    % Apply impulsive control if it's an impulse time
    if is_impulse_time_good(k+1)
        for i = 1:n
            impulse_sum = 0;
            for j = 1:n
                if A(i,j) > 0
                    % Note the sign here - matches the thesis equations
                    impulse_sum = impulse_sum + A(i,j) * (X_with_condition(j, k+1) - X_with_condition(i, k+1));
                end
            end
            X_with_condition(i, k+1) = X_with_condition(i, k+1) + impulsive_strength * impulse_sum;
        end
    end
end

% Simulation without necessary conditions satisfied
for k = 1:N-1
    % Calculate control input for each agent
    for i = 1:n
        u_i = 0;
        for j = 1:n
            if A(i,j) > 0
                diff_ij = X_without_condition(j, k) - X_without_condition(i, k);
                u_i = u_i + A(i,j) * (c1 * sgn_alpha(diff_ij) + ...
                                     c2 * sgn_gamma(diff_ij) + ...
                                     c3 * sign(diff_ij));
            end
        end
        
        % Update state with continuous dynamics
        X_without_condition(i, k+1) = X_without_condition(i, k) + dt * u_i;
    end
    
    % Apply impulsive control if it's an impulse time
    if is_impulse_time_bad(k+1)
        for i = 1:n
            impulse_sum = 0;
            for j = 1:n
                if A(i,j) > 0
                    % For the destabilizing case, we need to ensure the sign is correct
                    impulse_sum = impulse_sum + A(i,j) * (X_without_condition(j, k+1) - X_without_condition(i, k+1));
                end
            end
            % Apply the impulse with strength > 1
            X_without_condition(i, k+1) = X_without_condition(i, k+1) + impulsive_strength * impulse_sum;
        end
    end
end

% Plot results with necessary conditions satisfied
figure(1);
plot(t, X_with_condition, 'LineWidth', 2);
grid on;
xlabel('t (seconds)');
ylabel('State Values (Node No = Initial Value)');
title('LLFxTC in case of impulsive strength greater than one - with necessary conditions');
legend('State 1', 'State 2', 'State 3', 'State 4');
axis([0 T -3 7]);  % Set proper axis limits

% Find the approximate consensus value
consensus_value = mean(X_with_condition(:, end));
% Add text annotation for consensus point
hold on;
text(0.3, consensus_value-0.5, ['CONSENSUS POINT', newline, 'x = 0.5', newline, 'y = ', num2str(consensus_value, '%.2f')], 'FontWeight', 'bold', 'Color', 'red');
hold off;

% Plot results without necessary conditions satisfied
figure(2);
plot(t, X_without_condition, 'LineWidth', 2);
grid on;
xlabel('t (seconds)');
ylabel('State Values (Node No = Initial Value)');
title('LLFxTC in case of impulsive strength greater than one - without necessary conditions');
legend('State 1', 'State 2', 'State 3', 'State 4');
axis([0 T -40 10]);  % Set proper axis limits to match the diverging behavior

% Visualize the graph topology
figure(3);
G = graph(A);
p = plot(G, 'LineWidth', 2);
p.NodeLabel = {'1', '2', '3', '4'};
title('Graph Topology in Example 2');
