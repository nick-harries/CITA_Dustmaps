import sys
import os
import numpy as np
from scipy.optimize import minimize
import re



""" 
This preamble calls the functions needed to optimize the array, which are located in separate directories. This should be updated so that the directory paths are located in their own file.

The preamble creates the frequencies array, bounds list, calls the necessary functions, and reads the fits files into a large array.
The bounds are defined in the data_processing code and can be changed there. The only variable that might need to be changed on this code is the 'nside', which is the resolution that you degrade the arrays to.
"""



# Read the file that you have set up with your directory paths
tex_file_path = 'directories.tex'
with open(tex_file_path, 'r') as file:
    tex_content = file.read()
config_dir_match = re.search(r"Configuration files directory:\s*'(.+?)'", tex_content)
fits_dir_match = re.search(r"Fits files directory:\s*'(.+?)'", tex_content)
# Define the directories
if config_dir_match and fits_dir_match:
    config_directory = config_dir_match.group(1)
    fits_directory = fits_dir_match.group(1)
else:
    raise ValueError("Directory paths not found in the tex file.")
# Add the configuration directory to the system path
sys.path.append(config_directory)

# Import the necessary functions from the associated codes.
from data_processing import extract_info, extract_constants_and_parameters, stokes_arrays_in_MJy_sr, decrease_resolution, conversion_factors_Kcmb_to_MJy_sr, blackbody_function, stokes_reconstruction, Chi2, define_bounds, extract_number
from plot_functions import plot_recreated_values_and_fit, plot_optimized_parameters, plot_optimized_parameters_histograms, plot_minimized_Chi2

# List all FITS files in the directory and create 'files' array to be used to extract stokes arrays
fits_files = [file for file in os.listdir(fits_directory) if file.endswith('.fits')]
fits_files = sorted(fits_files, key = extract_number)
files = [os.path.join(fits_directory, file) for file in fits_files]


# Define the path to the constants and parameters file, then extract constants and parameters arrays
constants_and_parameters_file = os.path.join(config_directory, 'constants_and_parameters.tex')
parameters, constants = extract_constants_and_parameters(constants_and_parameters_file)

#Define array of bounds as are set in the data_processing.py function
bounds = define_bounds()

#Create array of frequency values associated to the fits files
frequencies = np.zeros(len(files))
for f in range(len(files)):
    frequencies[f] = np.float64(extract_info(files[f])[6])
frequencies = np.array([1.43e11, 2.17e11, 3.43e11, 5.45e11, 8.57e11])
print(frequencies)

#Define nside resolution to change stokes arrays to, then extract information from fits files, change units to MJy/sr, change resolution, and store array



nside = 128
print('nside: ', nside)
nest_type = extract_info(files[0])[-1]

print('creating arrays ...')
stokes_arrays_full_resolution = stokes_arrays_in_MJy_sr(files, constants, frequencies)
stokes_arrays_correct_units = decrease_resolution(stokes_arrays_full_resolution, nest_type, nside)
print('arrays created ...')


def execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies):
    """
    This function employs the scipy.optimize.minimize function to optimize the parameters for each pixel. This is done by finding the set of parameters that minimizes the returned value of the Chi^2 function that is defined above.
    *This function is currently set up to simply take a sample of the real data, to test its efficiency*
    *The starting and ending indices of this sample are also returned, to aid in plotting*
    Two empty arrays are created, stokes_arrays_reconstructed, which has shape (n , 3, l), where n is the number of fits files the code has read (equivalent to the number of frequencies being modeled), 3 is hard-coded because stokes I, Q, and U will be recreated for each frequency, l is the length of the arrays within the fits file, this is dependent on the chosen nside resolution.
    optimized_parameters_array has shape (l, 7), where l is the same as above, and 7 is hard coded for the 6 parameters being optimized (Temperature, Tau, alpha, Beta, Psi, p_frac) and the final index stores the final Chi2 value from the optimizer. Each optimized parameter for each pixel is stored in this array.
    As the code is currently set up, each of the scipy.optimize.minimize methods Nelder-Mead, BFGS, TNC, and Powell are attempted, and that which returns the smallest minimized Chi2 value has its optimized parameters stored. If the method returns an error, its output is ignored.
    The arrays_to_optimize list is created, storing which index of the stokes data array actually contains data rather than being full of nans (This will happen if the fits file does not have data for a stokes parameter), this information is passed to the Chi2 calculator to ensure Chi2 is calculated based off of the frequencies and parameters with real data.
    The tally_count is simply to see how many results were returned by which scipy.optimize.minimize method, this has no effect on the optimization.
    
    Within the for loop, i refers to the pixel within the Planck data that we are optimizing, j refers to the index of the empty array that is being filled with data, all j values can be replaced with i if starting from the 0th pixel or optimizing the entire array.
    Initial guesses for parameters is kept constant as successively updating seems to preferentially optimize stokes I data.

    Input: Tuple of constants, stokes_arrays_correct units, bounds, parameters, array of frequencies
    Output: Array of reconstructed stokes parameters (I, Q, U for each frequency input into the function), Array of optimized parameters and minimized Chi2 value corresponding to each optimized pixel, list of which index of the stokes data array contains real data (Used for plotting), starting index of sample(used for plotting), ending index of sample (used for plotting)
    """

    #Create list of which indices contain numerical data and should be used to calculate Chi2

    print('create coordinates list...')
    arrays_to_optimize = []
    for file_number in range(np.shape(stokes_arrays_correct_units)[0]):
        for parameter in range(3):
            if not np.isnan(stokes_arrays_correct_units[file_number, parameter, 0]):
                arrays_to_optimize.append([file_number, parameter])
    print(arrays_to_optimize)   
   

    length_of_sample = int(len(stokes_arrays_correct_units[0, 0, :]))
    starting_index = int(0)
    ending_index = int(starting_index + length_of_sample)

    stokes_arrays_reconstructed = np.zeros((len(frequencies), 3, length_of_sample))
    optimized_parameters_array = np.zeros((length_of_sample, 7))

    optimize_methods = ['BFGS', 'TNC', 'Powell']
    tally_count = [0, 0, 0, 0]

    print('begin optimization, starting index: %d, ending index: %d ...'%(starting_index, ending_index))
    j = 0
    for i in range(starting_index, ending_index): #Range of i is the range of indices of the real data that we will optimize parameters for

        optimized_params = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method='Nelder-Mead')#, bounds=bounds)
        minimized_Chi_2 = optimized_params.fun

        for idx, method in enumerate(optimize_methods):
            try:
                # Alternate optimization with different methods, see which returns smallest Chi2 value
                alternate_optimization = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method=method) #, bounds=bounds)
                if alternate_optimization.fun < minimized_Chi_2:
                    minimized_Chi_2 = alternate_optimization.fun
                    optimized_params = alternate_optimization
                    tally_count[idx + 1] += 1  #Update tally count for which method was used
            except: #Ignore error messages
                continue



        recreated_values = stokes_reconstruction(optimized_params.x, constants, frequencies) #Models emission based off of the optimized parameters, returns I, Q, U values in each input frequency
        # Store the emission reconstruction values in the empty array
        for freq in range(np.shape(recreated_values)[0]): #For each input frequency
            for param in range(np.shape(recreated_values)[1]): #For each parameter (I, Q, U)
                stokes_arrays_reconstructed[freq, param, j] = recreated_values[freq, param]
        
        for param in range(6): # Store each of the 6 optimized parameters, will be in the same order as the initial parameters array (T, Beta, Tau, Psi, Alpha, p_frac)
            optimized_parameters_array[j, param] = optimized_params.x[param]
        optimized_parameters_array[j, 6] = minimized_Chi_2 #Store the minimal Chi2 value for the pixel in the last column of the array

        print(j, '%.3e'%minimized_Chi_2)
        j+=1 #Update index number for storing in arrays

    tally_count[0] = length_of_sample - tally_count[1] - tally_count[2] - tally_count[3] #Number of times the initial method is used is the length of the array minus the number of times the other methods were used.
    print('tally count: ', tally_count)

    return stokes_arrays_reconstructed, optimized_parameters_array, arrays_to_optimize, starting_index, ending_index

stokes_arrays_reconstructed, optimized_parameters_array, arrays_to_optimize, starting_index, ending_index = execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies)
print(np.shape(stokes_arrays_reconstructed))
print(np.shape(optimized_parameters_array))
np.save('Stokes_Arrays_CMB_Removed_nside_%d'%nside, stokes_arrays_reconstructed)
np.save('Parameters_Arrays_CMB_Removed_nside_%d'%nside, optimized_parameters_array)

print('Plotting ...')
plot_recreated_values_and_fit(stokes_arrays_reconstructed, stokes_arrays_correct_units, arrays_to_optimize, starting_index, ending_index, frequencies)
#plot_optimized_parameters(optimized_parameters_array)
plot_minimized_Chi2(optimized_parameters_array)
#plot_optimized_parameters_histograms(optimized_parameters_array)
