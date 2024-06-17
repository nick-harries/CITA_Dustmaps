import sys
import os
import numpy as np
from scipy.optimize import minimize

# Define the directory containing configuration files and data_processing.py
config_directory = '/home/nicholas-harries/Desktop/Planck_copies/configuration_files'
fits_directory = '/home/nicholas-harries/Desktop/Planck_copies/fits_files'

# Add the configuration directory to the system path
sys.path.append(config_directory)

# Now you can import the functions from data_processing
from data_processing import extract_info, extract_constants_and_parameters, stokes_arrays_in_MJy_sr, decrease_resolution, conversion_factors_Kcmb_to_MJy_sr, blackbody_function, stokes_reconstruction, Chi2, define_bounds
from plot_functions import plot_recreated_values_and_fit, plot_optimized_parameters, plot_optimized_parameters_histograms, plot_raw_data_against_itself, plot_raw_data, plot_minimized_Chi2

# List all FITS files in the directory
fits_files = [file for file in os.listdir(fits_directory) if file.endswith('.fits')]

# Construct full paths to the FITS files
files = [os.path.join(fits_directory, file) for file in fits_files]


# Define the path to the constants and parameters file
constants_and_parameters_file = os.path.join(config_directory, 'constants_and_parameters.tex')

# Now you can use the full path to the constants and parameters file
files = (files[0], files[1])  # Use the full paths
print(files)


bounds = define_bounds()


# Extract frequencies using the correct file paths
frequencies = np.array([extract_info(files[0])[6], extract_info(files[1])[6]], dtype=np.float64)
parameters, constants = extract_constants_and_parameters(constants_and_parameters_file)
stokes_arrays_correct_units = stokes_arrays_in_MJy_sr(files, constants, frequencies)
stokes_arrays_correct_units = decrease_resolution(stokes_arrays_correct_units, 32)

def execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters):
    """
    This function employs the scipy.optimize.minimize function to optimize the parameters for each pixel. This is done by finding the set of parameters that minimizes the returned value of the Chi^2 function that is defined above.
    *This function is currently set up to simply take a sample of the real data, to test its efficiency*
    *The starting and ending indices of this sample are also returned, to aid in plotting*
    An empty array is initialized, this array has the shape (n, 13), where n is the number of pixels that will be optimized in a given run.
    Once the optimized parameters are found, the stokes emission is recreated, these values (I 217GHz, Q 217GHz, U 217GHz, I 353GHz, Q 353GHz, U 353GHz) are stored in the first 5 columns of the empty array.
    The next 6 columns of the empty array are filled with the optimized parameters for the specific pixel in the order of Temperature, Beta, Tau, Psi, Alpha, p_frac.
    The 13th column of the empty array is filled with the minimized Chi^2 value, this is useful for trying different optimization methods to see which works best.

    Within the for loop, i refers to the pixel within the Planck data that we are optimizing, j refers to the index of the empty array that is being filled with data, all j values can be replaced with i if starting from the 0th pixel or optimizing the entire array.
    Initial guesses for parameters is kept constant as successively updating seems to preferentially optimize stokes I data.

    Input: Tuple of constants, stokes_arrays_correct units, bounds, parameters
    Output: Array of shape (n,13), where n is the number of pixels that have been optimized. This array contains modelled emission based off of the optimized parameters in both frequencies, the optimized parameters corresponding to each pixel, and the optimized Chi^2 value for each pixel.
    """
    # Initialize array to store the data
    print('length of array: ', len(stokes_arrays_correct_units[0, 0, :]))
    length_of_array = len(stokes_arrays_correct_units[0, 0, :])
    length_of_sample = length_of_array // 6
    data = np.zeros((length_of_sample, 13))  # First number is the number of iterations, 13 columns: 6 columns for 6 recreated emission values, 6 columns for optimized parameters, 1 column for optimized Chi^2 value
    j = 0
    starting_index = 0 * length_of_sample
    ending_index = 1 * length_of_sample
    for i in range(starting_index, ending_index): #Range of i is the range of indices of the real data that we will optimize parameters for

        optimized_params = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, i), method='L-BFGS-B', bounds=bounds)
        minimized_Chi_2 = optimized_params.fun


        recreated_values = stokes_reconstruction(optimized_params.x, constants, frequencies) #Models emission based off of the optimized parameters, returns I, Q, U values in 217Ghz and 353Ghz

        # Store the emission reconstruction values in the empty array
        data[j, 0] = recreated_values[0, 0]  # Stokes I 217GHz
        data[j, 1] = recreated_values[0, 1]  # Stokes Q 217GHz
        data[j, 2] = recreated_values[0, 2]  # Stokes U 217GHz
        data[j, 3] = recreated_values[1, 0]  # Stokes I 353GHz
        data[j, 4] = recreated_values[1, 1]  # Stokes Q 353GHz
        data[j, 5] = recreated_values[1, 2]  # Stokes U 353GHz

        # Store the minimized parameters in the empty array
        data[j, 6] = optimized_params.x[0] # Temperature
        data[j, 7] = optimized_params.x[1] # Beta
        data[j, 8] = optimized_params.x[2] # Tau
        data[j, 9] = optimized_params.x[3] # Psi
        data[j, 10] = optimized_params.x[4] # Alpha
        data[j, 11] = optimized_params.x[5] #p_frac

        #Store the minimized output of the Chi^2 function
        data[j, 12] = minimized_Chi_2
        j+=1

    return data, starting_index, ending_index

data, starting_index, ending_index = execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters)

print(data)

plot_recreated_values_and_fit(data, stokes_arrays_correct_units, starting_index, ending_index)
plot_optimized_parameters(data)
plot_minimized_Chi2(data)
plot_optimized_parameters_histograms(data)