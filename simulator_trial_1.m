%% 2-Link Planar Manipulator with Nonlinear Tendon Extensions
clc; clear; close all;

%% Link Geometry and Mass Properties (Cuboids)
% Link 1 dimensions
L1 = 0.45; % [m]
W1 = 0.18; % [m]
H1 = 0.15; % [m]
m1 = 0.012;  % [kg]

floor = true;

% Link 2 dimensions
L2 = 0.25;  % [m]
W2 = 0.13; % [m]
H2 = 0.11; % [m]
m2 = 0.007;  % [kg]

% Gravity
g = -0.981; % [m/s^2]
% g = 0; % [m/s^2]

%% Moments of Inertia via Cuboid Formula
% Moment of inertia about the z-axis through CM: I_z(cm) = (1/12)*m*(h^2 + w^2)

I1_cm = (1/12)*m1*(H1^2 + W1^2);
d1 = L1/2; 
I1 = I1_cm + m1*d1^2; % parallel axis shift for link 1

I2_cm = (1/12)*m2*(H2^2 + W2^2);
% For link 2, we use I2_cm directly in the standard formula for a 2-link arm
I2 = I2_cm;

%% Tendon Parameters
% From the snippet, four tendons: h1, h2, h3, h4.
l1 = 0.3; 
l2 = 0.3;
l3 = 0.35;
l4 = 0.35;
a = 0.005; % example value
b = 0.003; % example value

R1 = 0.002;  
R2 = 0.002;

% Tendon forces
%f = [0.007355; 0.007355; 0.007355; 0.007355]; % [N]
% f = @(t) [0.007355*exp(-5*t); 0.007355*exp(-5*t); 0.007355*exp(-5*t); 0.007355*exp(-5*t)];
% % f = @(t) [0.005*exp(-5*t); 0.007*exp(-5*t); 0.007*exp(-5*t); 0.007*exp(-5*t)];
% f = @(t) [0.007*exp(-5*t); 0.000*exp(-5*t); 0.000*exp(-5*t); 0.007*exp(-5*t)];
% f = @(t) [0.7*exp(-5*t); 0.00*exp(-5*t); 0.00*exp(-5*t); 0.7*exp(-5*t)];
% f = @(t) [0.07; 0.007355; 0.007355; 0.07];
% f = @(t) [0; 0; 7; 8.5];
f = @(t) [0; 0; 100*exp(-4*t); 85*exp(-10*t^2)];


%% Kinematics: Nonlinear Tendon Extension Functions
% According to the snippet (for θ1 > 0):
% h1(θ) = l1 + 2√(a²+b²)*cos(tan^-1(a/b) + θ1/2) - 2b - R2 θ2
% h2(θ) = l2 - R1 θ1
% h3(θ) = l3 + R1 θ1
% h4(θ) = l4 + R1 θ1 + R2 θ2
% If θ1 < 0, relations may be reversed, but we focus on θ1 ≥ 0 scenario.

H = @(theta) [
    l1 + 2*sqrt(a^2+b^2)*cos(atan(a/b) + theta(1)/2) - 2*b - R2*theta(2);  % h1
    l2 - R1*theta(1);                                                     % h2
    l3 + R1*theta(1);                                                     % h3
    l4 + R1*theta(1) + R2*theta(2)                                        % h4
];

% Compute P(theta) = dH/dtheta:
% ∂h1/∂θ1 = -2*√(a²+b²)*sin(...) * (1/2) = -√(a²+b²)*sin(...)
% ∂h1/∂θ2 = -R2
% ∂h2/∂θ1 = -R1; ∂h2/∂θ2 = 0
% ∂h3/∂θ1 = R1;  ∂h3/∂θ2 = 0
% ∂h4/∂θ1 = R1;  ∂h4/∂θ2 = R2

dH_dtheta = @(theta) [
    -2*sqrt(a^2+b^2)*sin(atan(a/b) + theta(1)/2),  0;
    -R1,                                         0;
    R1,                                          -R2;
    R1,                                          R2
];

% P(theta) = (dH/dtheta)^T
P_theta = @(theta) dH_dtheta(theta)';

%% Dynamics (M, C, N)
% Define link centers of mass
lc1 = L1/2;
lc2 = L2/2;

% Standard parameters from 2-link arm dynamics
a_m = m1*lc1^2 + m2*(L1^2+lc2^2) + I1 + I2;
b_m = m2*L1*lc2;
d_m = m2*lc2^2 + I2;

M_theta = @(th) [a_m+2*b_m*cos(th(2)), d_m+b_m*cos(th(2));
                 d_m+b_m*cos(th(2)),   d_m];

C_theta = @(th, thd) [ -b_m*sin(th(2))*thd(2), -b_m*sin(th(2))*(thd(1)+thd(2));
                        b_m*sin(th(2))*thd(1), 0 ];

N_theta = @(th) [-(m1*g*lc1*cos(th(1)) + m2*g*L1*cos(th(1)))+m2*g*lc2*cos(th(1)-th(2));
                 -m2*g*lc2*cos(th(1)-th(2))];

% Equation of motion:
% M(theta)*ddtheta + C(theta,theta_dot)*theta_dot + N(theta) = P(theta)*f
% => ddtheta = M(theta)\( P(theta)*f - C(theta,theta_dot)*theta_dot - N(theta) )

theta_ddot_func = @(th, thd, t) M_theta(th)\(P_theta(th)*f(t) - C_theta(th, thd)*thd - N_theta(th));

%% Simulation Parameters
tmax =.5;  
dt = 0.001;
n = round(tmax/dt);
t = linspace(0, tmax, n);

% Initial conditions
theta0 = [deg2rad(-30); deg2rad(-90)]; 
theta_dot0 = [0; 0];

theta = zeros(2, n);
theta_dot = zeros(2, n);
theta_ddot = zeros(2, n);

theta(:,1) = theta0;
theta_dot(:,1) = theta_dot0;

%% Numerical Integration (Euler's Method)

% for i = 1:n-1
%     f_current = f(t(i)); % Compute time-varying force
%     theta_ddot(:, i) = M_theta(theta(:, i))\(P_theta(theta(:, i))*f_current ...
%                        - C_theta(theta(:, i), theta_dot(:, i))*theta_dot(:, i) ...
%                        - N_theta(theta(:, i)));
%     theta_dot(:, i+1) = theta_dot(:, i) + theta_ddot(:, i)*dt;
%     theta(:, i+1) = theta(:, i) + theta_dot(:, i)*dt;
% end

theta_ddot(:, end) = theta_ddot_func(theta(:, end), theta_dot(:, end), t(end));

% Desired limits in radians
if floor == true
    theta1_min = deg2rad(-30);
else
    theta1_min = deg2rad(-75);
end
theta1_max = deg2rad(10);
theta2_min = deg2rad(-90);
theta2_max = deg2rad(0);

for i = 1:n-1
    f_current = f(t(i));
    theta_ddot(:, i) = M_theta(theta(:, i))\(P_theta(theta(:, i))*f_current ...
                       - C_theta(theta(:, i), theta_dot(:, i))*theta_dot(:, i) ...
                       - N_theta(theta(:, i)));
    theta_dot(:, i+1) = theta_dot(:, i) + theta_ddot(:, i)*dt;
    theta_temp = theta(:, i) + theta_dot(:, i)*dt;

    % Clamp joint 1 angle
    if theta_temp(1) > theta1_max
        theta_temp(1) = theta1_max;
        theta_dot(1, i+1) = 0;
    elseif theta_temp(1) < theta1_min
        theta_temp(1) = theta1_min;
        theta_dot(1, i+1) = 0;
    end

    % Clamp joint 2 angle
    if theta_temp(2) > theta2_max
        theta_temp(2) = theta2_max;
        theta_dot(2, i+1) = 0;
    elseif theta_temp(2) < theta2_min
        theta_temp(2) = theta2_min;
        theta_dot(2, i+1) = 0;
    end

    theta(:, i+1) = theta_temp;
end


%% Plot Results
figure;
subplot(3,1,1);
plot(t, theta(1,:), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('\theta_1 [rad]');
title('Joint 1 Angle');
grid on;

subplot(3,1,2);
plot(t, theta(2,:), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('\theta_2 [rad]');
title('Joint 2 Angle');
grid on;

subplot(3,1,3);
plot(t, theta_dot(1,:), 'LineWidth', 2); hold on;
plot(t, theta_dot(2,:), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Joint Angular Velocities');
legend('\theta_1 dot', '\theta_2 dot');
grid on;

%% Animation of the Manipulator
figure('Name', '2-Link Manipulator Animation');
axis equal;
xlim([-1, 1]);
ylim([-1, 1]);
xlabel('X [m]');
ylabel('Y [m]');
grid on; hold on;

link1_plot = line([0,0],[0,0],'LineWidth',4,'Color','b');
link2_plot = line([0,0],[0,0],'LineWidth',4,'Color','r');

for i=1:10:n
    th1 = theta(1,i);
    th2 = theta(2,i);

    x1 = L1*cos(th1);
    y1 = L1*sin(th1);
    x2 = x1 + L2*cos(th1+th2);
    y2 = y1 + L2*sin(th1+th2);

    set(link1_plot, 'XData',[0,x1],'YData',[0,y1]);
    set(link2_plot, 'XData',[x1,x2],'YData',[y1,y2]);

    pause(.01);
    drawnow;
end

%% Save Animation as MP4 Video
videoFileName = '2LinkManipulatorSimulation.mp4'; % Desired video file name
videoWriter = VideoWriter(videoFileName, 'MPEG-4'); % Use MPEG-4 profile for MP4
videoWriter.FrameRate = 20; % Adjust frame rate as needed
open(videoWriter);

% Animation Loop
figure('Name', '2-Link Manipulator Animation');
axis equal;
xlim([-1, 1]);
ylim([-1, 1]);
xlabel('X [m]');
ylabel('Y [m]');
grid on; hold on;

link1_plot = line([0,0],[0,0],'LineWidth',4,'Color','b');
link2_plot = line([0,0],[0,0],'LineWidth',4,'Color','r');

for i=1:10:n
    th1 = theta(1,i);
    th2 = theta(2,i);

    x1 = L1*cos(th1);
    y1 = L1*sin(th1);
    x2 = x1 + L2*cos(th1+th2);
    y2 = y1 + L2*sin(th1+th2);

    set(link1_plot, 'XData',[0,x1],'YData',[0,y1]);
    set(link2_plot, 'XData',[x1,x2],'YData',[y1,y2]);

    % Capture the current frame
    frame = getframe(gcf);
    writeVideo(videoWriter, frame);

    pause(.001);
    drawnow;
end

% Close the video writer
close(videoWriter);

disp(['Simulation saved as ', videoFileName]);
