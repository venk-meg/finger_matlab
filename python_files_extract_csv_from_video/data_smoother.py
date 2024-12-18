import pandas as pd
import numpy as np

# Load the data
filename = "folddown/corrected_data.csv"  # Replace with your file name
data = pd.read_csv(filename)

def smooth_data(data, window_size=5):
    """
    Smooth the data using a moving average filter.
    
    Parameters:
        data (pd.DataFrame): Input data with landmark columns.
        window_size (int): The size of the moving average window.
        
    Returns:
        pd.DataFrame: Smoothed data.
    """
    smoothed_data = data.copy()
    for joint in ['knuckle', 'PIP', 'DIP', 'tip']:
        for axis in ['x', 'y', 'z']:
            column = f"{joint}_{axis}"
            # Apply moving average filter
            smoothed_data[column] = smoothed_data[column].rolling(window=window_size, center=True).mean()
    return smoothed_data

# Apply smoothing
window_size = 5  # Adjust the window size as needed
smoothed_data = smooth_data(data, window_size)

# Handle NaN values introduced by the moving average at the edges
smoothed_data.fillna(method='bfill', inplace=True)  # Backward fill
smoothed_data.fillna(method='ffill', inplace=True)  # Forward fill

# Save the smoothed data
output_filename = "folddown/smoothed_data.csv"
smoothed_data.to_csv(output_filename, index=False)

print(f"Smoothed data saved to {output_filename}")
