% Load the CSV data
data = readtable('folddown/planar_data_best_fit.csv');
%video saved as
vidname = 'folddown/planar_data_best_fit_finger_simulation_adjusted.mp4'
%figure saved as
figname = 'folddown/planar_data_best_fit_finger_simulation_adjusted.fig'

% Extract the coordinates
knuckle = [data.knuckle_x, data.knuckle_y, data.knuckle_z];
PIP = [data.PIP_x, data.PIP_y, data.PIP_z];
DIP = [data.DIP_x, data.DIP_y, data.DIP_z];
tip = [data.tip_x, data.tip_y, data.tip_z];

% Set up the figure
fig = figure; % Store the figure handle
hold on;
grid on;
axis equal; % Equal aspect ratio for all axes
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('4-Link Finger Simulation in 3D');

% Compute the maximum range to ensure all points fit symmetrically around the origin
all_points = [knuckle; PIP; DIP; tip]; % Combine all points
max_extent = max(max(abs(all_points), [], 1)); % Maximum distance from the origin in any direction

% Set symmetric axis limits around the origin
xlim([-max_extent, max_extent]);
ylim([-max_extent, max_extent]);
zlim([-max_extent, max_extent]);

% Plot the X, Y, and Z axes for clarity
line([-max_extent, max_extent], [0, 0], [0, 0], 'Color', 'k', 'LineStyle', '--'); % X-axis
line([0, 0], [-max_extent, max_extent], [0, 0], 'Color', 'k', 'LineStyle', '--'); % Y-axis
line([0, 0], [0, 0], [-max_extent, max_extent], 'Color', 'k', 'LineStyle', '--'); % Z-axis

axis manual; % Fix the axis scale
view(3); % 3D view

% Initialize MP4 video writer
videoFile = VideoWriter(vidname, 'MPEG-4');
videoFile.FrameRate = 30;
open(videoFile);

% Animate the finger motion
for i = 1:height(data)
    % Extract positions for the current frame
    points = [knuckle(i, :); PIP(i, :); DIP(i, :); tip(i, :)];
    
    % Clear previous plot
    cla;
    
    % Plot the links in 3D
    plot3(points(:,1), points(:,2), points(:,3), '-o', ...
          'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    
    % Plot the X, Y, and Z axes for clarity
    line([-max_extent, max_extent], [0, 0], [0, 0], 'Color', 'k', 'LineStyle', '--'); % X-axis
    line([0, 0], [-max_extent, max_extent], [0, 0], 'Color', 'k', 'LineStyle', '--'); % Y-axis
    line([0, 0], [0, 0], [-max_extent, max_extent], 'Color', 'k', 'LineStyle', '--'); % Z-axis

    % Set axis labels and grid
    xlabel('X (meters)'); ylabel('Y (meters)'); zlabel('Z (meters)');
    grid on;
    
    % Capture the current frame
    frame = getframe(gcf);
    writeVideo(videoFile, frame);
    
    % Pause to simulate video playback speed (30 FPS)
    pause(1 / 30); % 1 / frame rate
end

% Close the video writer
close(videoFile);

% Save the figure as an interactive .fig file
savefig(fig, figname);

disp('Animation saved');
disp('Interactive 3D figure saved');
