import pandas as pd
import numpy as np

# Load the Illuminator data
file_path = '/Users/defneatac/Documents/UCL_MSc/MSc Disso/Illuminator_hem/Scenarios/pv_data_Rotterdam_NL-15min.txt'
illuminator_data = pd.read_csv(file_path, delimiter=',')

# Extract the timestamps from the Illuminator data using the correct column name
time_column = 'Time'  # Adjust based on the actual column name in the file
illuminator_timestamps = pd.to_datetime(illuminator_data[time_column])
 
# Create a DataFrame for the HEM data
hem_data = {
    "SimulationTime": {
        "start": 0,
        "end": 24,
        "step": 1
    },
    "ExternalConditions": {
        "air_temperatures": [1.0] * 24,
        "wind_speeds": [3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
        "ground_temperatures": [8.0, 8.7, 9.4, 10.1, 10.8, 10.5, 11.0, 12.7, 8.0, 8.7, 9.4, 10.1, 10.8, 10.5, 11.0, 12.7, 8.0, 8.7, 9.4, 10.1, 10.8, 10.5, 11.0, 12.7],
        "diffuse_horizontal_radiation": [0, 0, 0, 0, 0, 0, 0, 0, 52, 85, 95, 88, 80, 58, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "direct_beam_radiation": [0, 0, 0, 0, 0, 0, 0, 0, 169, 284, 328, 397, 354, 205, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "solar_reflectivity_of_ground": [0.2, 0.22, 0.2, 0.22, 0.247, 0.24, 0.26, 0.233, 0.287, 0.233, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.213, 0.227, 0.22, 0.24],
    }
}
 
# Create a DataFrame for HEM data
hem_df = pd.DataFrame({
    'Time': pd.date_range(start='2012-01-01 00:00:00', periods=24, freq='H'),
    'air_temperatures': hem_data["ExternalConditions"]["air_temperatures"],
    'wind_speeds': hem_data["ExternalConditions"]["wind_speeds"],
    'wind_directions': hem_data["ExternalConditions"]["ground_temperatures"],  # Assuming ground_temp as a placeholder for wind_directions
    'diffuse_horizontal_radiation': hem_data["ExternalConditions"]["diffuse_horizontal_radiation"],
    'direct_beam_radiation': hem_data["ExternalConditions"]["direct_beam_radiation"],
    'solar_reflectivity_of_ground': hem_data["ExternalConditions"]["solar_reflectivity_of_ground"]
})
 
# print(hem_df)
 
# Extend the HEM DataFrame to match the time range of the Illuminator data
extended_hem_df = hem_df.set_index('Time').reindex(pd.date_range(start=illuminator_timestamps.min(), end=illuminator_timestamps.max(), freq='15min'))
 
# Interpolate the extended HEM data to fill in missing values
extended_hem_df = extended_hem_df.interpolate(method='time')
 
# Debug: Print the HEM data after interpolation
print("Extended HEM data after interpolation:\n", extended_hem_df.head(10))
 
# Write the interpolated data to a new text file with the required format

# Write the interpolated data to a new text file with the required format
output_file_path = '/Users/defneatac/Documents/UCL_MSc/MSc Disso/Illuminator_hem/Scenarios/interpolated_external_conditions.txt'
extended_hem_df.head(100).to_csv(output_file_path, index_label='Time,', header=True)

print(f"Interpolated data written to {'/Users/defneatac/Documents/UCL_MSc/MSc Disso/Illuminator_hem/Scenarios'}")

