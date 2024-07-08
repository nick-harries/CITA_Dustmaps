import sys
import os
import numpy as np
from scipy.optimize import minimize
import re
import time

# Read the file that you have set up with your directory paths
txt_file_path = 'my_directories.txt'
with open(txt_file_path, 'r') as file:
    txt_content = file.read() #Search within the directories.txt file
config_dir_match = re.search(r"Configuration files directory:\s*'(.+?)'", txt_content) #Find configuration file path
fits_dir_match = re.search(r"Fits files directory:\s*'(.+?)'", txt_content) #Find Fits file directory path

#Troubleshoot to make sure that the directory paths exist before executing rest of code
if config_dir_match and fits_dir_match:
    config_directory = config_dir_match.group(1)
    fits_directory = fits_dir_match.group(1)
else: #If the directory paths do not exist, print error message and stop code
    raise ValueError("Directory paths not found in the tex file.")
    exit()
# Add the configuration directory to the system path
sys.path.append(config_directory)

# Import the necessary functions from the associated codes.
from data_processing import extract_info, extract_constants_and_parameters, stokes_arrays_in_MJy_sr, decrease_resolution, conversion_factors_Kcmb_to_MJy_sr, blackbody_function, stokes_reconstruction, Chi2, define_bounds, extract_number
from plot_functions import plot_recreated_values_and_fit, plot_optimized_parameters, plot_optimized_parameters_histograms, plot_minimized_Chi2


""" 
This preamble calls the functions needed to optimize the array, which are located in separate directories. This should be updated so that the directory paths are located in their own file.

The preamble creates the frequencies array, bounds list, calls the necessary functions, and reads the fits files into a large array.
The bounds are defined in the data_processing code and can be changed there. The only variable that might need to be changed on this code is the 'nside', which is the resolution that you degrade the arrays to.
"""

# List all FITS files in the directory and create 'files' array to be used to extract stokes arrays
fits_files = [file for file in os.listdir(fits_directory) if file.endswith('.fits')]
fits_files = sorted(fits_files, key = extract_number) #Frequency maps from PLA have frequencies in the name, sort through files and open them in order of increasing frequency
#The line above is necessary because CMB-removed maps do not have frequency values within the metadata
files = [os.path.join(fits_directory, file) for file in fits_files]


# Define the path to the constants and parameters file, then extract constants and parameters arrays
constants_and_parameters_file = os.path.join(config_directory, 'constants_and_parameters.txt')
parameters, constants = extract_constants_and_parameters(constants_and_parameters_file)

#Define array of bounds as are set in the data_processing.py function
bounds = define_bounds()

#Create array of frequency values associated to the fits files
frequencies = np.zeros(len(files))
for f in range(len(files)):
    frequencies[f] = extract_info(files[f])[6] #Sort through the fits files and append frequency values extracted from metadata
if any(np.isnan(frequencies)): #CMB-removed maps have no frequency information in their metadata, this condition assumes frequencies associated with files in that case, should be modified to fit your condition
    print('Frequencies not in metadata, using hard-coded values')
    frequencies = np.array([143., 217., 343., 545., 857.])
print('Files attributes...')
for freq in range(frequencies.shape[0]):
    print(f'{frequencies[freq]} GHz, {extract_info(files[freq])[-1]}')

#Define nside resolution to change stokes arrays to, then extract information from fits files, change units to MJy/sr, change resolution, and store the modified array
nside = 128
print('Converting to nside: ', nside)
nest_type = extract_info(files[0])[-2]

print('creating arrays ...')
stokes_arrays_full_resolution = stokes_arrays_in_MJy_sr(files, constants, frequencies)
stokes_arrays_correct_units = decrease_resolution(stokes_arrays_full_resolution, nest_type, nside)
print('arrays created ...')


def execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies):

    """
    This function employs the scipy.optimize.minimize function to optimize the parameters for each pixel. This is done by finding the set of parameters that minimizes the returned-
    -value of the Chi^2 function that is defined above.
    *This function is currently set up to simply take a sample of the real data, to test its efficiency*
    *The starting and ending indices of this sample are also returned, to aid in plotting*
    Two empty arrays are created, stokes_arrays_reconstructed, which has shape (n , 3, l), where n is the number of fits files the code has read (equivalent to the number of frequencies being modeled),-
    - 3 represents stokes I, Q, and U, which will be recreated for each frequency, l is the length of the arrays within the fits file, this is dependent on the chosen nside resolution.
    optimized_parameters_array has shape (l, 7), where l is the same as above, and is for the 6 parameters being optimized (Temperature, Tau, alpha, Beta, Psi, p_frac)-
    - and the final index stores the final Chi2 value from the optimizer. Each optimized parameter for each pixel is stored in this array.
    As the code is currently set up, each of the scipy.optimize.minimize methods Nelder-Mead, BFGS, and Powell are attempted, and that which returns the smallest minimized Chi2 -
    - value has its optimized parameters stored. If the method returns an error, its output is ignored.
    The arrays_to_optimize list is created, storing which index of the stokes data array actually contains data rather than being full of zeros (for stokes parameters) or ones (for covariances) -
    - (This will happen if the fits file does not have data for a stokes parameter), this information is passed to the Chi2 calculator to ensure Chi2 is calculated based -
    - off of the frequencies and parameters with real data.
    The tally_count is simply to see how many results were returned by which scipy.optimize.minimize method, this has no effect on the optimization.
    
    Within the for loop, i refers to the pixel within the Planck data that we are optimizing, j refers to the index of the empty array that is being filled with data, -
    - all j values can be replaced with i if starting from the 0th pixel or optimizing the entire array.
    Initial guesses for parameters is kept constant as successively updating seems to preferentially optimize stokes I data.

    Input: 
        Constants - This is extracted from the constants and parameters text file by the extract_constants_and_parameters function in data_processing.py. It is a tuple and is used in the blackbody function used to recreate stokes parameters
        stokes_arrays_correct_units - This is the multidimensional array of data from the imported fits files. the Shape is (a, 6, c), where a is the number of files that were read, 6 is for the 3 stokes parameters and their associated covariances. If any data for a parameter or covariance is not present in the fits file its column is set as either ones or zeros. c is the length of an array of whatever nside value the map was deresolved to.
        bounds - These are the limits that the parameters will be optimized within. The bounds are set in the define_bounds function in data_processing.py
        parameters - This is the tuple of fiducial parameter values that is read from the constants and parameters text file by the extract_constants_and_parameters function in the data_processing.py code. These values serve as the initial guess for each pixel.
        frequencies - This is the array of frequency values associated with each of the fits files that was read. It is used to recreate the stokes emission.
    Output: 
        stokes_arrays_reconstructed - This is the multidimensional array of modeled emission. Its shape is (a, 3, b), where a is the number of frequencies across which the optimizer is running, and b is the length of the sample of the data that is being optimized. If the entire map is being optimized, b will be equal to the length of an array of the given nside. 3 is for the 3 stokes parameters that were modeled, I, Q, and U (in that order).
        optimized_parameters_array - This is the multidimensional array of optimized emission parameters for each pixel. Its shape is (c, 7), where c is the number of pixels (length of an array of the given nside), and 7 is for the 6 emission parameters, (T, Beta, Tau, Psi, Alpha, p_frac), with the 7th index being the final Chi2 value for the pixel.
        arrays_to_optimize - This is an array of shape (a, b) containing the coordinates within stokes_arrays_correct_units array that contain actual data (i.e, Planck 545GHz maps do not contain polarization data). It is used for the plotting functions.
        starting_index - This is an integer value that tells the first index of the healpy arrays that was optimized in the current sample. This will be set to 0 for an optimization of the entire map. This is used in plotting functions to line up the data indices with the model indices.
        ending_index - This is an integer value that tells the last index of the healpy arrays that was optimized in the current sample. This will be set to the length of an array of the given nside if the entire map is being plotted. This is used in plotting functions to specify the range of data values that are plotted against the modeled emission.
    """

    #Define arrays containing the elements of the model that will be optimized and stored
    stokes_parameters = np.array(['I', 'Q', 'U'])
    parameters_list = np.array(['Temperature', 'Beta', 'Tau', 'Psi', 'Alpha', 'p_frac', 'Chi2'])

    print('Optimization bounds: ') 
    for index, param in enumerate(parameters_list[:-1]):
        print(f'{bounds[index][0]} <= {param} <= {bounds[index][1]}')


    #Create array of which indices contain numerical data and should be used to calculate Chi2
    print('Finding array coordinates with data...')
    arrays_to_optimize = []
    for file_number in range(stokes_arrays_correct_units.shape[0]): #For each frequency map
        for parameter in range(stokes_parameters.shape[0]): #For each of the three stokes parameters
            if not np.all(stokes_arrays_correct_units[file_number, parameter, :] == 0): #If the row of emission data is not all zeros (Checking that the given frequency contains the specified parameter)
                arrays_to_optimize.append([file_number, parameter]) #Append the coordinates within stokes_arrays_correct_units to the list, only if the coordinates contain nonzero data
    arrays_to_optimize = np.array(arrays_to_optimize) #Convert to array

    #The coordinates to be used for the Chi2 value are now stored. This print statement tells you which parameters exist for which frequencies
    for coord in arrays_to_optimize:
        a, b = coord
        print(f"{frequencies[a]} GHz", stokes_parameters[b])
   

    length_of_sample = 100#Choose how many indices to optimize in this sample
    starting_index = int(3e4)#Which index of the array to start at
    ending_index = starting_index + length_of_sample #Which index to end at (automatically calculated from length_of_sample and starting_index variables

    #The following arrays are returned at the end of the function, containing the optimized values
    stokes_arrays_reconstructed = np.zeros((frequencies.shape[0], stokes_parameters.shape[0], length_of_sample)) #initialize empty array to store stokes array models
    optimized_parameters_array = np.zeros((length_of_sample, parameters_list.shape[0])) #Initialize empty array to store optimized parameters for each pixel of each optimized array

    optimize_methods_with_bounds = np.array(['L-BFGS-B', 'trust-constr']) #scipy.optimize.minimize methods, all are used to find best possible parameter combination
    optimize_methods_no_bounds = np.array(['Nelder-Mead', 'BFGS', 'TNC'])
    tally_count =np.zeros(optimize_methods_no_bounds.shape[0]+1)#Used to see which optimization methods are most efficient, stores the number of times each method is used

    print(f'begin optimization, starting index: {starting_index}, ending index: {ending_index} ...')
    j = 0 #j indexes within the arrays of stored data, i indexes within the data arrays
    start_time = time.time() #Timer started right before optimizing loop is called
    for i in range(starting_index, ending_index): #Range of i is the range of indices of the real data that we will optimize parameters for

        optimization_results = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method='Powell')# bounds=bounds)
        optimized_params = optimization_results.x #Store the tuple of optimized parameters
        minimized_Chi_2 = optimization_results.fun #Store the minimized Chi2 value

        for idx, method in enumerate(optimize_methods_no_bounds): #Cycle through optimizer methods, if new method returns smaller Chi2 value, store its parameters, ignore output in case of ValueErrors
            try:
                # Alternate optimization with different methods, see which returns smallest Chi2 value
                alternate_optimization = minimize(Chi2, parameters, args=(frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i), method=method)# bounds=bounds)
                if alternate_optimization.fun < minimized_Chi_2: #If the returned Chi2 value is smaller than the previously stored version, keep it and its associated parameters
                    minimized_Chi_2 = alternate_optimization.fun
                    optimized_params = alternate_optimization.x
                    tally_count[idx + 1] += 1  #Update tally count for which method was used
            except: #Ignore error messages, somewhat common with optimizers
                continue


        # Recreate values using the optimized parameters
        recreated_values = stokes_reconstruction(optimized_params, constants, frequencies)  # Models emission based on the optimized parameters, returns I, Q, U values in each input frequency

        #j is the pixel index. Store the recreated stokes parameters for each of the frequencies being optimized at the current pixel index
        stokes_arrays_reconstructed[:, :, j] = recreated_values

        # Store each of the 6 optimized parameters
        optimized_parameters_array[j, :-1] = optimized_params

        # Store the minimal Chi2 value for the pixel in the last column of the array
        optimized_parameters_array[j, -1] = minimized_Chi_2

        print(f'j: {j}, Chi2: {minimized_Chi_2}') #Print how far into the length of the sample you are, as well as the final Chi2 value for the jth pixel. Delete this line for efficiency when running a long sample
        j+=1 #Update index number for storing in arrays


    end_time = time.time() #Store the time directly after the entire sample is done being optimized
    elapsed_time = end_time - start_time #The amount of time the optimizer took to run is the difference between the time values before and after it was initiated
    time_for_whole_array = (stokes_arrays_correct_units.shape[2] / length_of_sample) * elapsed_time #Extrapolate how long it would have taken to optimize the entire map rather than just a sample

    # Convert time_for_whole_array to hours, minutes, and seconds
    hours = int(time_for_whole_array // 3600)
    minutes = int((time_for_whole_array % 3600) // 60)
    seconds = int(time_for_whole_array % 60)

    #Print indicators of the success of the run, Chi2 values, time taken
    tally_count[0] = length_of_sample - np.sum(tally_count) #The first element of the tally_count array is for the 'Nelder-Mead' method and is not updated. This line calculates the number of times it was used
    print('Tally count: ', tally_count)
    print(f'Time taken for sample: {elapsed_time:.2f} seconds')
    print(f'Time estimated for whole array: {hours}h {minutes}m {seconds}s')
    print(f'Mean Chi2: {np.mean(optimized_parameters_array[:, -1])}')

    return stokes_arrays_reconstructed, optimized_parameters_array, arrays_to_optimize, starting_index, ending_index

stokes_arrays_reconstructed, optimized_parameters_array, arrays_to_optimize, starting_index, ending_index = execute_Chi2_optimization(constants, stokes_arrays_correct_units, bounds, parameters, frequencies)
print(np.shape(stokes_arrays_reconstructed))
print(np.shape(optimized_parameters_array))
#np.save('Stokes_Arrays_CMB_Removed_nside_%d'%nside, stokes_arrays_reconstructed)
#np.save('Parameters_Arrays_CMB_Removed_nside_%d'%nside, optimized_parameters_array)

print('Plotting ...')
plot_recreated_values_and_fit(stokes_arrays_reconstructed, stokes_arrays_correct_units, arrays_to_optimize, starting_index, ending_index, frequencies)
plot_optimized_parameters(optimized_parameters_array)
plot_minimized_Chi2(optimized_parameters_array)
#plot_optimized_parameters_histograms(optimized_parameters_array)
