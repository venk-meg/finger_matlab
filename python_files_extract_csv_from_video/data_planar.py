# import pandas as pd
# import numpy as np
# from sklearn.decomposition import PCA

# # Load the data
# filename = "folddown/smoothed_data.csv"  # Replace with your file name
# data = pd.read_csv(filename)

# def fit_and_project_to_plane(data):
#     """
#     Fit a plane to the points and project all points onto it.

#     Parameters:
#         data (pd.DataFrame): Input data with landmark columns.

#     Returns:
#         pd.DataFrame: Data projected onto the best-fit plane.
#     """
#     projected_data = data.copy()

#     for frame in range(len(data)):
#         # Combine all points for this frame
#         points = []
#         for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
#             points.append(data.loc[frame, [f"{joint}_x", f"{joint}_y", f"{joint}_z"]].to_numpy())
#         points = np.array(points)

#         # Fit a plane using PCA
#         pca = PCA(n_components=2)
#         pca.fit(points)
#         plane_normal = np.cross(pca.components_[0], pca.components_[1])  # Normal vector to the plane

#         # Project points onto the plane
#         for i, joint in enumerate(['knuckle', 'PIP', 'DIP', 'tip']):
#             point = points[i]
#             distance_to_plane = np.dot(point, plane_normal) / np.linalg.norm(plane_normal)
#             projected_point = point - distance_to_plane * plane_normal
#             projected_data.loc[frame, [f"{joint}_x", f"{joint}_y", f"{joint}_z"]] = projected_point

#     return projected_data

# # Project data onto the best-fit plane
# planar_data = fit_and_project_to_plane(data)

# # Save the planar data
# output_filename = "folddown/planar_data.csv"
# planar_data.to_csv(output_filename, index=False)

# print(f"Planar data saved to {output_filename}")


import pandas as pd
import numpy as np

# Load the data
filename = "folddown/smoothed_data.csv"  # Replace with your file name
data = pd.read_csv(filename)


def force_planarity(data, target_plane='xy'):
    """
    Force all points to lie on a single plane (e.g., X-Y plane or best-fit plane).
    
    Parameters:
        data (pd.DataFrame): Input data with landmark columns.
        target_plane (str): 'xy' for X-Y plane (Z=0), or 'best-fit' for data-driven plane.
        
    Returns:
        pd.DataFrame: Adjusted data with points forced onto the target plane.
    """
    adjusted_data = data.copy()

    if target_plane == 'xy':
        # Project all points onto the X-Y plane (Z = 0)
        for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
            adjusted_data[f"{joint}_z"] = 0

    elif target_plane == 'best-fit':
        # Combine all points from all frames
        all_points = []
        for frame in range(len(data)):
            for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
                all_points.append(data.loc[
                    frame,
                    [f"{joint}_x", f"{joint}_y", f"{joint}_z"]].to_numpy())
        all_points = np.array(all_points)

        # Compute the best-fit plane using PCA
        mean_point = np.mean(all_points, axis=0)
        centered_points = all_points - mean_point
        _, _, vh = np.linalg.svd(centered_points)
        plane_normal = vh[-1]  # Normal vector of the best-fit plane
        print(plane_normal)

        # Project each point onto the best-fit plane
        for frame in range(len(data)):
            for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
                point = data.loc[
                    frame,
                    [f"{joint}_x", f"{joint}_y", f"{joint}_z"]].to_numpy()
                distance_to_plane = np.dot(point - mean_point, plane_normal)
                projected_point = point - distance_to_plane * plane_normal
                adjusted_data.loc[frame,
                                  [f"{joint}_x", f"{joint}_y", f"{joint}_z"
                                   ]] = projected_point

    return adjusted_data

# Force the data onto the X-Y plane or best-fit plane
planar_data_xy = force_planarity(data, target_plane='xy')  # Force onto X-Y plane
planar_data_best_fit = force_planarity(data, target_plane='best-fit')  # Force onto best-fit plane

# Save the adjusted data
output_filename_xy = "folddown/planar_data_xy.csv"
output_filename_best_fit = "folddown/planar_data_best_fit.csv"
planar_data_xy.to_csv(output_filename_xy, index=False)
planar_data_best_fit.to_csv(output_filename_best_fit, index=False)

print(f"Planar data saved to {output_filename_xy} and {output_filename_best_fit}")
