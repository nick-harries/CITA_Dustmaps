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
from plot_functions import plot_recreated_values_and_fit, plot_optimized_parameters, plot_optimized_parameters_histograms, plot_minimized_Chi2

# List all FITS files in the directory
fits_files = [file for file in os.listdir(fits_directory) if file.endswith('.fits')]

# Construct full paths to the FITS files
files = [os.path.join(fits_directory, file) for file in fits_files]

# Define the path to the constants and parameters file
constants_and_parameters_file = os.path.join(config_directory, 'constants_and_parameters.tex')

bounds = define_bounds()

frequencies = np.zeros(len(files))

for f in range(len(files)):
    frequencies[f] = np.float64(extract_info(files[f])[6])

print(frequencies)
# Extract frequencies using the correct file paths
parameters, constants = extract_constants_and_parameters(constants_and_parameters_file)
stokes_arrays_correct_units = stokes_arrays_in_MJy_sr(files, constants, frequencies)
stokes_arrays_correct_units = decrease_resolution(stokes_arrays_correct_units, 32)

def execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies):
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
    
    arrays_to_optimize = []
    for file_number in range(np.shape(stokes_arrays_correct_units)[0]):
        for parameter in range(np.shape(stokes_arrays_correct_units)[1]):
            if stokes_arrays_correct_units[file_number, parameter, 0] is not np.nan:
                arrays_to_optimize.append([file_number, parameter])

    print(arrays_to_optimize)
   
   
    length_of_sample = int(100)
    starting_index = int(1.5e3)
    ending_index = int(starting_index + length_of_sample)

    stokes_arrays_reconstructed = np.zeros((len(frequencies), 3, length_of_sample))
    optimized_parameters_array = np.zeros((length_of_sample, 7))

    optimize_methods = ['Nelder-Mead', 'TNC', 'Powell', 'Trust-Constr']

    j = 0
    for i in range(starting_index, ending_index): #Range of i is the range of indices of the real data that we will optimize parameters for

        optimized_params = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method='BFGS')#, bounds=bounds)
        minimized_Chi_2 = optimized_params.fun

        for method in optimize_methods:
            try:
                alternate_optimization = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method=optimize_methods[method])#, bounds=bounds)
                if alternate_optimization.fun < minimized_Chi_2:
                    minimized_Chi_2 = alternate_optimization.fun
                    optimized_params = alternate_optimization
            except:
                continue


        recreated_values = stokes_reconstruction(optimized_params.x, constants, frequencies) #Models emission based off of the optimized parameters, returns I, Q, U values in 217Ghz and 353Ghz
        # Store the emission reconstruction values in the empty array
        for freq in range(np.shape(recreated_values)[0]):
            for param in range(np.shape(recreated_values)[1]):
                stokes_arrays_reconstructed[freq, param, j] = recreated_values[freq, param]
        
        for param in range(6):
            optimized_parameters_array[j, param] = optimized_params.x[param]
        optimized_parameters_array[j, 6] = minimized_Chi_2

        print(j)
        j+=1

    return stokes_arrays_reconstructed, optimized_parameters_array, starting_index, ending_index

stokes_arrays_reconstructed, optimized_parameters_array, starting_index, ending_index = execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies)


plot_recreated_values_and_fit(stokes_arrays_reconstructed, stokes_arrays_correct_units, starting_index, ending_index, frequencies)
#plot_optimized_parameters(optimized_parameters_array)
plot_minimized_Chi2(optimized_parameters_array)
#plot_optimized_parameters_histograms(optimized_parameters_array)