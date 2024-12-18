import pandas as pd
import numpy as np

# Load the data
filename = "folddown/folddown_iphone_index_finger_landmarks_test_20241130_205754_adjusted.csv"  # Replace with your file name
data = pd.read_csv(filename)

def compute_rotation_matrix(v):
    """
    Compute the rotation matrix to align vector `v` with the X-Y plane.
    """
    v = v / np.linalg.norm(v)  # Normalize the vector
    # Project `v` onto the Z=0 plane
    projection = np.array([v[0], v[1], 0])
    projection = projection / np.linalg.norm(projection)  # Normalize projection
    
    # Compute the rotation axis (cross product of `v` and projection)
    axis = np.cross(v, projection)
    axis = axis / np.linalg.norm(axis)  # Normalize axis

    # Compute the angle between `v` and its projection
    angle = np.arccos(np.dot(v, projection))
    
    # Use Rodrigues' rotation formula to compute the rotation matrix
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    
    return R

def rotate_points(points, R):
    """
    Rotate a set of points using rotation matrix `R`.
    """
    rotated_points = np.dot(points, R.T)  # Apply the rotation matrix
    return rotated_points

def scale_points(data, scaling_factor):
    """
    Scale the points in the data by the given scaling factor.
    """
    for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
        for axis in ['x', 'y', 'z']:
            data[f"{joint}_{axis}"] *= scaling_factor
    return data

def apply_180_degree_rotation(data):
    """
    Apply a 180-degree rotation around the x-axis to the entire dataset.
    """
    R_180_x = np.array([
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ])
    for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
        points = data[[f'{joint}_x', f'{joint}_y', f'{joint}_z']].to_numpy()
        rotated_points = rotate_points(points, R_180_x)
        data[[f'{joint}_x', f'{joint}_y', f'{joint}_z']] = rotated_points
    return data

# Step 1: Find the frame where tip is approximately at [-0.0377, -0.0025, 0.0751]
target_tip = np.array([-0.0377, -0.0025, 0.0751])
closest_frame_idx = np.argmin(
    np.linalg.norm(data[['tip_x', 'tip_y', 'tip_z']].to_numpy() - target_tip, axis=1)
)

# Step 2: Scale the data
# Compute the scaling factor using the closest frame (knuckle to PIP distance = 45mm)
knuckle = data.loc[closest_frame_idx, ['knuckle_x', 'knuckle_y', 'knuckle_z']].to_numpy()
PIP = data.loc[closest_frame_idx, ['PIP_x', 'PIP_y', 'PIP_z']].to_numpy()
actual_distance = np.linalg.norm(PIP - knuckle)
scaling_factor = 45 / actual_distance

# Apply scaling to all points
data = scale_points(data, scaling_factor)

# Step 3: Rotate the data to align with the X-Y plane
# Compute the rotation matrix for the closest frame
tip = data.loc[closest_frame_idx, ['tip_x', 'tip_y', 'tip_z']].to_numpy()
vector = tip - knuckle
rotation_matrix = compute_rotation_matrix(vector)

# Apply rotation to all points
for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
    points = data[[f'{joint}_x', f'{joint}_y', f'{joint}_z']].to_numpy()
    rotated_points = rotate_points(points, rotation_matrix)
    data[[f'{joint}_x', f'{joint}_y', f'{joint}_z']] = rotated_points

# Step 4: Apply 180-degree rotation around the X-axis
data = apply_180_degree_rotation(data)

# Save the scaled and rotated data
output_filename = "folddown/corrected_data.csv"
data.to_csv(output_filename, index=False)

print(f"Scaled, rotated, and 180-degree inverted data saved to {output_filename}")
