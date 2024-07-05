from astropy.io import fits
import healpy as hp
import numpy as np
import re

def extract_info(fits_file):

    """
    This function calls the specified fits file and goes through its HDU, storing information useful for plotting and other calculations.
    In particular, it calls and returns any available stokes parameters(I,Q,U), and associated covariances.
    It also stores the Ordering type for healpix plotting as well as the frequency of the file. If the frequency is not in units of GHz, an error message is printed and the code is stopped.
    Soft-Coded so it will not crash if any information is not available.

    Input is the string name of the fits file, assuming it is in the same directory as the code
    Output is stokes i, stokes q, stokes u, covariance of i, covariance of q, covariance of u, frequency of file, condition about ordering type
    """
    
    with fits.open(fits_file) as hdul:
        if len(hdul) > 1:
            if 'ORDERING' in hdul[1].header:
                if hdul[1].header['ORDERING'] == 'NESTED':
                    nest_type = True
                else:
                    nest_type = False
            else:
                print('nest type of %s unknown...'%fits_file)
                exit

            if hdul[1].data is not None:
                data = hdul[1].data
               
                if 'I_STOKES' in data.dtype.names:
                    i_stokes_data = data['I_STOKES']
                    i_stokes_data[i_stokes_data < 0] = 0 #This line resets physically impossible negative values to 0
                elif 'INTENSITY' in data.dtype.names:
                    i_stokes_data = data['INTENSITY'].flatten()
                    i_stokes_data[i_stokes_data < 0] = 0

                else:
                    i_stokes_data = None

                if 'Q_STOKES' in data.dtype.names:
                    q_stokes_data = data['Q_STOKES']
                elif 'Q-POLARISATION' in data.dtype.names:
                    q_stokes_data = data['Q-POLARISATION'].flatten()
                else:
                    #q_stokes_data = None
                    q_stokes_data = np.zeros_like(i_stokes_data)
                if 'U_STOKES' in data.dtype.names:
                    u_stokes_data = data['U_STOKES']
                elif 'U-POLARISATION' in data.dtype.names:
                    u_stokes_data = data['U-POLARISATION'].flatten()
                else:
                    #u_stokes_data = None
                    u_stokes_data = np.zeros_like(i_stokes_data)
                if 'II_COV' in data.dtype.names:
                    ii_cov_data = data['II_COV']
                else:
                    #ii_cov_data = None
                    ii_cov_data = np.ones_like(i_stokes_data)
                if 'UU_COV' in data.dtype.names:
                    uu_cov_data = data['UU_COV']
                else:
                    #uu_cov_data = None
                    uu_cov_data= np.ones_like(i_stokes_data)
                if 'QQ_COV' in data.dtype.names:
                    qq_cov_data = data['QQ_COV']
                else:
                    #qq_cov_data = None 
                    qq_cov_data = np.ones_like(i_stokes_data)
            if 'FREQ' in hdul[1].header:
                freq = int(hdul[1].header['FREQ'])
                if 'UNITFREQ' in hdul[1].header:
                    if hdul[1].header['UNITFREQ'] == 'GHz':
                        freq = freq * 10 ** 9
                    else:
                        print('Frequency unit error, not GHz. Fatal error')
                        exit
            else:
                freq = None

    return i_stokes_data, q_stokes_data, u_stokes_data, ii_cov_data, qq_cov_data, uu_cov_data, freq, nest_type
    
def extract_constants_and_parameters(const_and_params_file):
    """ 
This function calls the text file containing the constants and parameters, and returns them as tuples.

Input: Name of text file containing constants and parameters as string
Returns: Tuple of parameters, Tuple of constants
    """

# Initialize empty lists for parameters and constants
    parameters = []
    constants = []

    # Open the file and read line by line
    with open(const_and_params_file, 'r') as file:
        for line in file:
            # Skip comment lines and empty lines
            if line.strip().startswith('#') or not line.strip():
                continue
        
            # Process the parameters lines
            if 'Parameter:' in line:
                # Extract value from the line
                value = line.split('Parameter:')[1].split('#')[0].strip()
                parameters.append(float(value))
        
            # Process the constants lines
            if 'Constant:' in line:
                # Extract value from the line
                value = line.split('Constant:')[1].split('#')[0].strip()
                constants.append(float(value))

    # Convert lists to tuples
    parameters = tuple(parameters)
    constants = tuple(constants)
    return parameters, constants

def stokes_arrays_in_MJy_sr(files, constants, frequencies):

    """
    This function calls the fits files and constants, and converts the stokes parameters and covariance arrays from kcmb/kcmb^2 to (MJy/sr)/(MJy/sr)^2. 
    First, it uses the constants and frequencies of the files to compute conversion factors between the units, then returns the converted values.
    For time-saving sake, the returned values are all in a single array of shape (n, 6, l), where n is the number of files being read, 6 corresponds to the 3 stokes parameters (I, Q, U) and their associated covariances, and l is the length of the arrays within the file.
    This is soft-coded so that if a fits file does not contain an array for a specific parameter, the parameter is filled in with NAN indices.

    Input: list/array/tuple of fits file names as strings, tuple of constants, array of frequencies of the files being input(externally created)
    Output: Array of shape (n, 6, l) containing stokes parameters and associated covariances at n different frequencies, all converted from units of Kcmb to MJy/sr, columns of NAN in the place of any parameter with missing data
    """


    """ 
    The layout of the data array is:
    [i, 0, :]. [i, 3, :] = stokes I array for frequency i, stokes I covariance array for frequency i
    [i, 1, :], [i, 4, :] = stokes Q array for frequency i, stokes Q covariance array for frequency i
    [i, 2, :], [i, 5, :] = stokes U array for frequency i, stokes U covariance array for frequency i
    """

    stokes_parameters = ['I', 'Q', 'U', 'II', 'QQ', 'UU']
    files_length = len(files)
    array_length = len(extract_info(files[0])[0])
    num_arrays = 6
    extracted_data_array = np.zeros((files_length, 6, array_length))

    print('calculating Kcmb to MJy/sr conversion factors ...')
    conversion_factors = conversion_factors_Kcmb_to_MJy_sr(constants, frequencies)
    print(f'conversion factors calculated as: {conversion_factors}')

    for i in range(files_length):
        for j in range(3):
            extracted_data_array[i, j, :] = extract_info(files[i])[j] * conversion_factors[i]
        for j in range(3,6):
            extracted_data_array[i, j, :] = extract_info(files[i])[j] * conversion_factors[i] ** 2
    print('Array converted')
    #print('seeing which parameters exist and converting them...')


    #for i in range(files_length):  # For each file
    #    info = extract_info(files[i]) # Retrieve once per file
    #
    #    for j in range(3):  # For each stokes parameter (I, Q, U)
    #        if info[j] is not None:
    #            extracted_data_array[i, j, :] = info[j] * conversion_factors[i]
    #        else:
    #            extracted_data_array[i, j, :] = np.nan
    
    #    for k in range(3, 6):  # For each stokes parameter covariance (II, QQ, UU)
    #        if info[k] is not None:
    #            extracted_data_array[i, k, :] = info[k] * conversion_factors[i] ** 2
    #        else:
    #            extracted_data_array[i, k, :] = np.ones(array_length)


    return extracted_data_array

def decrease_resolution_original(extracted_data_array, nest_type, new_nside):
    """ 
    This function serves to decrease the resolution of the map extracted from the fits file. Less resolved maps have shorter arrays and can therefore be optimized quicker, for troubleshooting.
    The new_nside argument specifies the nside resolution that the outputted map will have.
    This function reads the first element of each column of the inputted array to see if it contains data or NAN (Which means the associated file/frequency does not have data for that parameter), only columns with numerical data will have resolution changed, empty columns remain empty.

    Input: Multidimensional array extracted from 'stokes_arrays_in_MJy_sr' function, desired nside value for that array
    Output: Multidimensional array of shape (n, 6, l), where n is the number of frequencies/files being used, 6 columns for the 3 stokes parameters and their associated covariances, l is the length of an array of nside=new_nside.
    """
    if nest_type == True:
        order='NESTED'
    else:
        order='RING'


    #Create blank array of shape (n, 6, l)
    array_length = hp.nside2npix(new_nside)
    extracted_data_array_new_nside = np.empty((len(extracted_data_array[:, 0, 0]), len(extracted_data_array[0, :, 0]), array_length))
    for frequency in range(extracted_data_array.shape[0]):
        for parameter in range(extracted_data_array.shape[1]):
            extracted_data_array_new_nside[frequency, parameter, :] = hp.ud_grade(extracted_data_array[frequency, parameter, :], new_nside, order_in = order, order_out = 'NESTED')
        
    return extracted_data_array_new_nside

def decrease_resolution(extracted_data_array, nest_type, new_nside):
    """ 
    This function serves to decrease the resolution of the map extracted from the fits file. Less resolved maps have shorter arrays and can therefore be optimized quicker, for troubleshooting.
    The new_nside argument specifies the nside resolution that the outputted map will have.
    This function reads the first element of each column of the inputted array to see if it contains data or NAN (Which means the associated file/frequency does not have data for that parameter), only columns with numerical data will have resolution changed, empty columns remain empty.

    Input: Multidimensional array extracted from 'stokes_arrays_in_MJy_sr' function, desired nside value for that array
    Output: Multidimensional array of shape (n, 6, l), where n is the number of frequencies/files being used, 6 columns for the 3 stokes parameters and their associated covariances, l is the length of an array of nside=new_nside.
    """
    order = 'NESTED' if nest_type else 'RING'

    # Create blank array of shape (n, 6, l)
    array_length = hp.nside2npix(new_nside)
    num_frequencies = extracted_data_array.shape[0]
    num_parameters = extracted_data_array.shape[1]
    
    extracted_data_array_new_nside = np.empty((num_frequencies, num_parameters, array_length))
    
    for frequency in range(num_frequencies):
        for parameter in range(num_parameters):
            data_slice = extracted_data_array[frequency, parameter, :]
            if not np.all(data_slice == 0):
                extracted_data_array_new_nside[frequency, parameter, :] = hp.ud_grade(data_slice, new_nside, order_in=order, order_out='NESTED')
            else:
                extracted_data_array_new_nside[frequency, parameter, :] = 0

    return extracted_data_array_new_nside

def conversion_factors_Kcmb_to_MJy_sr(constants, frequencies):

    """ 
    This function computes the conversion factor between units of kcmb and MJy/sr as a function of frequency. The tuple of constants is called for this, and the last value, freq_0=353e9GHz is explicitely ignored.
    The equation used to convert between frequencies is the temperature derivative of the Planck function, it is a function of frequency.

    Inputs: Array of frequency values. These should be previously defined before calling this function. Tuple of constants.
    Returns: Array of conversion factors corresponding to the elements within the array of frequencies.
    """
    c, h, k, T_cmb, *ignore = constants

#    T_cmb = 2.7255

    watt_m2_hz_sr_to_MJy_sr_conversion_factor = 10 ** 20 #The derivative of the Planck function converts from Kcmb to Watt / (m^2 Hz sr), this factor converts from Watt / (m^2 Hz sr) to MJy / sr

    x = (h * frequencies) / (k * T_cmb)

    kcmb_t0_watt_m2_hz_sr_factor = ( (2 * h * frequencies ** 3 )/ (c ** 2 * T_cmb) ) * ( (x * np.exp(x)) / (np.exp(x) - 1) ** 2)


    return kcmb_t0_watt_m2_hz_sr_factor * watt_m2_hz_sr_to_MJy_sr_conversion_factor

def blackbody_function(frequencies, constants, Temperature):
    """
This function calculates spectral radiance from Planck's law, as a function of Temperature and frequency. Temperature is a parameter that will be optimized in the Chi^2 function. Frequencies are an array that is previously defined.
The last two values of the 'constants' tuple, T_cmb and freq_0 are not needed and are explicitely ignored.
The spectral radiance is converted to correct MJy/sr units via the conversion factor 1(MJy/sr) = 10^20 (W/m^2 Hz sr)

Input: Array of frequencies, tuple of constants, Temperature parameter
Returns: Array of spectral radiance in units of MJy/sr at the same temperature value, elements of the array correspond to the elements of the frequencies array that is being called.
    """

    c, h, k_b, *ignore = constants #Speed of light(m/s), Planck's constant(J*s), Boltzmann's constant (J/K)
    #print(k_b)
    
    Watt_per_m2_Hz_sr_to_MJy_sr = 10 ** 20 #Planck's function calculates spectral radiance in units of Watt / (m^2 Hz sr), this factor converts from Watt / (m^2 Hz sr) to MJy / sr

    return Watt_per_m2_Hz_sr_to_MJy_sr * (2 * h * frequencies ** 3) / (c ** 2 * (np.exp((h * frequencies) / (k_b * Temperature)) - 1))
    
def stokes_reconstruction(parameters, constants, frequencies):

    """
    This function employs a distance-independant version of equations 2, 3, and 4 from Solaeche et al. 2018. 
    It models dust emission in stokes I, Q, and U parameters as a function of Temperature, spectral index (beta), optical depth (tau), frequency, and the parameters described below:
    p_frac is the polarization fraction, it applies only to the Q and U outputs.
    psi is an orientation term regarding the emitting dust and the magnetic field, it applies only to the Q and U ouputs
    alpha is a geometrical term regarding the direction of the magnetic field with respect to the line of sight, it applies only to the Q and U outputs
    k is arbitrarily defined as k=3, following Solaeche et al. 2018 and further, Fauvet et al. 2011

    Inputs: array of parameters, tuple of physical constants, array of frequencies
    Outputs: Recreated stokes emission in a single array of shape(2,3), where the first column corresponds to the columns of the input frequencies array and the second column corresponds to stokes I, Q, and U, respectively.
    """

    Temperature, beta, tau, psi, alpha, p_frac = parameters
    freq_0, k = constants[-1], constants[-2]

    stokes_i =  blackbody_function(frequencies, constants, Temperature) * tau * (frequencies / freq_0) ** beta
    stokes_q = stokes_i * np.cos(2 * psi) * (np.sin(alpha)) ** k * p_frac
    stokes_u = stokes_i * np.sin(2 * psi) * (np.sin(alpha)) ** k * p_frac

    return np.stack((stokes_i, stokes_q, stokes_u), axis=-1)

def Chi2(parameters, frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize, i):
    
    
    #recreated_values = stokes_reconstruction(parameters, constants, frequencies)
    """
    This function quantifies the fit of the modelled emission vs the emission from the data. The recreated emission is calculated in the stokes_reconstruction function.
    Within the 'execute_Chi2_optimization' function of the 'parameter_optimizer.py' code, a list called 'arrays_to_optimize' is created, this list contains the coordinates within the 'stokes_arrays_correct_units' -
    - array that contain numerical data. Then, the pixels within these columns are explicitely used to calculate the Chi2 value. Columns can be full of NAN indices if the fits file does not contain data for that parameter, this is meant to ignore those columns.
    Any stokes parameter column that contains data is assumed to contain an associated covariance column that is 3 indices farther down in the stokes_array_correct_units array. This shift is fundamental to how the arrays are extracted amd stored from the fits files so should not need to be changed despite its 'hard-coding'.
    The Chi2 value is calculated based off of the fit of the emission model with an initial set of parameters compared to the data from the fits files, the parameters are then optimized to find the smallest possible Chi2 value, thus the best fitting parameters.

    Inputs: Parameters array(will be optimized), frequencies, constants, stokes_arrays_correct_units, arrays_to_optimize list, i (i is the index of the loop that this function must be run through, it is defined in the 'execute_Chi2_optimization function)
    Outputs: Sum of all Chi^2 values for each stokes parameters across the array of frequencies
    """

    recreated_values = stokes_reconstruction(parameters, constants, frequencies)
    #print(arrays_to_optimize.shape)
    # Extract coordinates from arrays_to_optimize
    freq, param = arrays_to_optimize[:, 0], arrays_to_optimize[:, 1]

    
    # Calculate the Chi2 values
    #Chi2 = np.sum(np.abs(recreated_values[freq, param] - stokes_arrays_correct_units[freq, param, i]) ** 2 / stokes_arrays_correct_units[freq, param + 3, i])
    difference_squared = (recreated_values[freq, param] - stokes_arrays_correct_units[freq, param, i]) ** 2
    chi2_values = difference_squared / stokes_arrays_correct_units[freq, param + 3, i]
    Chi2 = np.sum(chi2_values)


    return Chi2

def define_bounds():
    """ 
    This is where the bounds that set the limits of optimization will be set.
    """
    return [(-np.infty, np.infty), #T
          (-np.infty, np.infty), #Beta
          (-np.infty, np.infty), #Tau
          (-np.infty, np.infty), #Psi(radians)
          (-np.infty, np.infty), #Alpha
          (-np.infty, np.infty)] #p_frac

def extract_number(filename):
    """
    Frequency map files  downloaded from the Planck Legacy Archive tend to have their associated frequency mentioned in the name, i.e, HFI_CompMap_Foregrounds-commander-143_R3.00.fits and HFI_SkyMap_143_2048_R3.01_full.fits are the 143GHz maps. This function finds the frequency from within the name of the file and sorts the fits files from smallest to highest frequency. This is implemented in case the frequency of the file is not in the metadata and thus the frequencies array must be hard-coded. This assures that the maps are in the same order as the hard-coded array.

    This function reads the name of the fits file, then stores the first sequence of numbers. It assumes that this number (143 in the examples above) is the frequency associated with the map. If a sequence is found, it is returned as an integer through the match.group(0) command, if no sequence is found within the title of the file name, a 0 is returned.
    """

    match = re.search(r'(\d+)', filename)
    return int(match.group(0)) if match else 0
