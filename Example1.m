% "Fixed Time Consensus Analysis of Multi Agent Systems via Impulsive Control"

clear all;
close all;

% Define the adjacency matrix A and Laplacian matrix L
A = [0 0.5 0 0 0;
     0.5 0 1.5 1 0;
     0 1.5 0 1 0;
     0 1 1 0 1;
     0 0 0 1 0];

L = [0.5 -0.5 0 0 0;
    -0.5 3 -1.5 -1 0;
     0 -1.5 2.5 -1 0;
     0 -1 -1 3 -1;
     0 0 0 -1 1];

% Parameters for the control protocol
c1 = 2;
c2 = 2;
c3 = 2;
alpha = 3/5;
gamma = 11/7;

% Initial states of the agents for Example 1
x0 = [5; -3; 6; -4; 1];

% Simulation parameters
T = 1.2; % Total simulation time
dt = 0.001; % Time step
t = 0:dt:T;
N = length(t);
n = length(x0); % Number of agents

% Function to calculate sgn^alpha(x)
sgn_alpha = @(x) sign(x).*abs(x).^alpha;

% Function to calculate sgn^gamma(x)
sgn_gamma = @(x) sign(x).*abs(x).^gamma;

% Array to store the states without impulsive control
X_without_impulse = zeros(n, N);
X_without_impulse(:, 1) = x0;

% Array to store the states with impulsive control
X_with_impulse = zeros(n, N);
X_with_impulse(:, 1) = x0;

% Impulsive parameters
impulsive_interval = 0.05; % Impulsive interval
impulsive_times = 0:impulsive_interval:T;
impulsive_strength = 0.75; % Impulsive strength (< 1)

% A more precise way to check impulsive times
is_impulse_time = zeros(1, N);
for i = 1:length(impulsive_times)
    [~, idx] = min(abs(t - impulsive_times(i)));
    is_impulse_time(idx) = 1;
end

% Simulation without impulsive control
for k = 1:N-1
    % Calculate control input for each agent
    for i = 1:n
        u_i = 0;
        for j = 1:n
            if A(i,j) > 0
                diff_ij = X_without_impulse(j, k) - X_without_impulse(i, k);
                u_i = u_i + A(i,j) * (c1 * sgn_alpha(diff_ij) + ...
                                     c2 * sgn_gamma(diff_ij) + ...
                                     c3 * sign(diff_ij));
            end
        end
        
        % Update state
        X_without_impulse(i, k+1) = X_without_impulse(i, k) + dt * u_i;
    end
end

% Simulation with impulsive control
for k = 1:N-1
    % Calculate control input for each agent
    for i = 1:n
        u_i = 0;
        for j = 1:n
            if A(i,j) > 0
                diff_ij = X_with_impulse(j, k) - X_with_impulse(i, k);
                u_i = u_i + A(i,j) * (c1 * sgn_alpha(diff_ij) + ...
                                     c2 * sgn_gamma(diff_ij) + ...
                                     c3 * sign(diff_ij));
            end
        end
        
        % Update state with continuous dynamics
        X_with_impulse(i, k+1) = X_with_impulse(i, k) + dt * u_i;
    end
    
    % Apply impulsive control if it's an impulse time
    if is_impulse_time(k+1)
        for i = 1:n
            impulse_sum = 0;
            for j = 1:n
                if A(i,j) > 0
                    impulse_sum = impulse_sum + A(i,j) * (X_with_impulse(j, k+1) - X_with_impulse(i, k+1));
                end
            end
            X_with_impulse(i, k+1) = X_with_impulse(i, k+1) + impulsive_strength * impulse_sum;
        end
    end
end

% Plot results without impulsive control
figure(1);
plot(t, X_without_impulse, 'LineWidth', 2);
grid on;
xlabel('t (seconds)');
ylabel('State Values (Node No = Initial Value)');
title('State Values vs t - Without Impulsive Control');
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5');
axis([0 T -5 7]);  % Set proper axis limits to see full plot

% Find the approximate consensus value
consensus_value = mean(X_without_impulse(:, end));
% Add text annotation for consensus point
hold on;
text(0.6, consensus_value, ['CONSENSUS POINT', newline, 'x = ', num2str(T), newline, 'y = ', num2str(consensus_value)], 'FontWeight', 'bold', 'Color', 'red');
hold off;

% Plot results with impulsive control
figure(2);
plot(t, X_with_impulse, 'LineWidth', 2);
grid on;
xlabel('t (seconds)');
ylabel('State Values (Node No = Initial Value)');
title('State Values vs t - With Impulsive Control (strength < 1)');
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5');
axis([0 T -2 6]);  % Set proper axis limits to see full plot

% Find the approximate consensus value
consensus_value_impulse = mean(X_with_impulse(:, end));
% Add text annotation for consensus point
hold on;
text(0.4, 3, ['CONSENSUS POINT', newline, 'x = 0.35', newline, 'y = ', num2str(consensus_value_impulse)], 'FontWeight', 'bold', 'Color', 'red');
hold off;

% Visualize the graph topology
figure(3);
G = graph(A);
p = plot(G, 'LineWidth', 2);
p.NodeLabel = {'1', '2', '3', '4', '5'};
p.EdgeLabel = G.Edges.Weight;
title('Graph Topology in Example 1');
