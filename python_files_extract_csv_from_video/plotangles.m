%% Load the CSV file
filename = 'data\smoothed_data.csv';
data = readtable(filename);

%% Extract positions
knuckle = [data.knuckle_x, data.knuckle_y, data.knuckle_z];
PIP = [data.PIP_x, data.PIP_y, data.PIP_z];
DIP = [data.DIP_x, data.DIP_y, data.DIP_z];
tip = [data.tip_x, data.tip_y, data.tip_z];

%% Time vector
time = (1:height(data))';

%% Calculate joint angles
% Angle at KNUCKLE (vector between KNUCKLE and PIP vs. horizontal axis)
k_to_p = PIP - knuckle;
k_to_d = DIP - PIP;
knuckle_angle = atan2d(vecnorm(cross(k_to_p, k_to_d, 2), 2, 2), dot(k_to_p, k_to_d, 2));

% Angle at PIP (vector between PIP and DIP vs. DIP to TIP)
p_to_d = DIP - PIP;
d_to_t = tip - DIP;
PIP_angle = atan2d(vecnorm(cross(p_to_d, d_to_t, 2), 2, 2), dot(p_to_d, d_to_t, 2));

%% Calculate angular velocities
knuckle_angle_velocity = [0; diff(knuckle_angle)]; % Numerical derivative
PIP_angle_velocity = [0; diff(PIP_angle)]; % Numerical derivative

%% Plot results
figure;

% Subplot 1: KNUCKLE angle vs time
subplot(3, 1, 1);
plot(time, knuckle_angle, 'LineWidth', 2);
xlabel('Time [frames]');
ylabel('KNUCKLE Angle [deg]');
title('KNUCKLE Angle vs Time');
grid on;

% Subplot 2: PIP angle vs time
subplot(3, 1, 2);
plot(time, PIP_angle, 'LineWidth', 2);
xlabel('Time [frames]');
ylabel('PIP Angle [deg]');
title('PIP Angle vs Time');
grid on;

% Subplot 3: Joint angular velocities
subplot(3, 1, 3);
plot(time, knuckle_angle_velocity, 'LineWidth', 2, 'DisplayName', 'KNUCKLE Velocity');
hold on;
plot(time, PIP_angle_velocity, 'LineWidth', 2, 'DisplayName', 'PIP Velocity');
legend('Location', 'best');
xlabel('Time [frames]');
ylabel('Angular Velocity [deg/frame]');
title('Joint Angular Velocities');
grid on;

%% Save the figure
saveas(gcf, 'joint_angles_and_velocities.png');
