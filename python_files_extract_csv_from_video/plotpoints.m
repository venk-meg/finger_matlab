% Load the CSV data
data = readtable('folddown/planar_data_best_fit.csv');
% Save figure name
figname_all_points_colored = 'folddown/planar_data_best_fit_all_points_colored_lines.fig';

% Extract the coordinates
knuckle = [data.knuckle_x, data.knuckle_y, data.knuckle_z];
PIP = [data.PIP_x, data.PIP_y, data.PIP_z];
DIP = [data.DIP_x, data.DIP_y, data.DIP_z];
tip = [data.tip_x, data.tip_y, data.tip_z];

% Combine all points for axis limits computation
all_points = [knuckle; PIP; DIP; tip];

% Set up the figure
fig_all_colored = figure; % Store the figure handle
hold on;
grid on;
axis equal; % Equal aspect ratio for all axes
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('3D Scatter Plot of All Finger Points with Connections');

% Compute the maximum range to ensure all points fit symmetrically around the origin
max_extent_all = max(max(abs(all_points), [], 1)); % Maximum distance from the origin in any direction

% Set symmetric axis limits around the origin
xlim([-max_extent_all, max_extent_all]);
ylim([-max_extent_all, max_extent_all]);
zlim([-max_extent_all, max_extent_all]);

% Plot the X, Y, and Z axes for clarity
line([-max_extent_all, max_extent_all], [0, 0], [0, 0], 'Color', 'k', 'LineStyle', '--'); % X-axis
line([0, 0], [-max_extent_all, max_extent_all], [0, 0], 'Color', 'k', 'LineStyle', '--'); % Y-axis
line([0, 0], [0, 0], [-max_extent_all, max_extent_all], 'Color', 'k', 'LineStyle', '--'); % Z-axis

% Plot each joint in a different color and connect with lines
for i = 1:height(data)
    % Extract the coordinates for the current frame
    points = [knuckle(i, :); PIP(i, :); DIP(i, :); tip(i, :)];
    
    % Plot the points
    scatter3(points(1,1), points(1,2), points(1,3), 50, 'r', 'filled', 'DisplayName', 'Knuckle');
    scatter3(points(2,1), points(2,2), points(2,3), 50, 'g', 'filled', 'DisplayName', 'PIP');
    scatter3(points(3,1), points(3,2), points(3,3), 50, 'b', 'filled', 'DisplayName', 'DIP');
    scatter3(points(4,1), points(4,2), points(4,3), 50, 'm', 'filled', 'DisplayName', 'Tip');
    
    % Connect the points with lines
    plot3(points(:,1), points(:,2), points(:,3), 'k-', 'LineWidth', 1.5);
end

% Add legend to differentiate joints
legend('Tip', 'DIP', 'PIP', 'Knuckle', 'Location', 'bestoutside');

% Set axis labels and grid
xlabel('X (meters)'); ylabel('Y (meters)'); zlabel('Z (meters)');
grid on;

% Save the figure as an interactive .fig file
savefig(fig_all_colored, figname_all_points_colored);

disp('3D Scatter Plot with connections saved as an interactive figure.');
