% Time settings
dt = 0.01; % Time step (s)
T = 5; % Total simulation time (s)
time = 0:dt:T; % Time vector

% Preallocate arrays for joint angles and velocities
theta = zeros(3, length(time));
theta_dot = zeros(3, length(time));

% Define robot parameters
L1 = 1; % Length of Link 1 (m)
L2 = 1; % Length of Link 2 (m)
L3 = 1; % Length of Link 3 (m)

% Initial conditions
theta0 = [0; 0; 0]; % Initial joint angles [rad]
theta_dot0 = [0; 0; 0]; % Initial angular velocities [rad/s]

% Define joint accelerations
theta_ddot = @(t) [0.1; 0.2; 0.3]; % Constant angular accelerations [rad/s^2]

% Initialize joint angles and velocities
theta(:, 1) = theta0;
theta_dot(:, 1) = theta_dot0;

% Integrate using Euler's method
for i = 2:length(time)
    % Update joint velocities
    theta_dot(:, i) = theta_dot(:, i-1) + theta_ddot(time(i)) * dt;

    % Update joint angles
    theta(:, i) = theta(:, i-1) + theta_dot(:, i) * dt;
end

% Preallocate arrays for end-effector positions
x = zeros(1, length(time));
y = zeros(1, length(time));

% Compute end-effector positions
for i = 1:length(time)
    x(i) = L1 * cos(theta(1, i)) + L2 * cos(theta(1, i) + theta(2, i)) + ...
           L3 * cos(theta(1, i) + theta(2, i) + theta(3, i));
    y(i) = L1 * sin(theta(1, i)) + L2 * sin(theta(1, i) + theta(2, i)) + ...
           L3 * sin(theta(1, i) + theta(2, i) + theta(3, i));
end

% Plot the end-effector trajectory
figure;
plot(x, y, 'r-', 'LineWidth', 2);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('End-Effector Trajectory');
grid on;


% Animate the robot motion
figure;
for i = 1:10:length(time)
    % Joint positions
    x1 = L1 * cos(theta(1, i));
    y1 = L1 * sin(theta(1, i));
    x2 = x1 + L2 * cos(theta(1, i) + theta(2, i));
    y2 = y1 + L2 * sin(theta(1, i) + theta(2, i));
    x3 = x2 + L3 * cos(theta(1, i) + theta(2, i) + theta(3, i));
    y3 = y2 + L3 * sin(theta(1, i) + theta(2, i) + theta(3, i));
    
    % Plot the robot
    plot([0, x1, x2, x3], [0, y1, y2, y3], 'bo-', 'LineWidth', 2);
    hold on;
    plot(x3, y3, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % End-effector
    hold off;
    axis equal;
    axis([-3, 3, -3, 3]);
    title('3-Link Planar Robot Motion');
    pause(0.05);
end

