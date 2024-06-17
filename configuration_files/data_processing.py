from astropy.io import fits
import healpy as hp
import numpy as np

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

            if hdul[1].data is not None:
                data = hdul[1].data

                if 'I_STOKES' in data.dtype.names:
                    epsilon = 1e-10
                    i_stokes_data = data['I_STOKES']
                    i_stokes_data[i_stokes_data < 0] = 0 #This line resets physically impossible negative values to 0
                else:
                    i_stokes_data = None

                if 'Q_STOKES' in data.dtype.names:
                    q_stokes_data = data['Q_STOKES']
                else:
                    q_stokes_data = None
                if 'U_STOKES' in data.dtype.names:
                    u_stokes_data = data['U_STOKES']
                else:
                    u_stokes_data = None

                if 'II_COV' in data.dtype.names:
                    ii_cov_data = data['II_COV']
                else:
                    ii_cov_data = None
                if 'UU_COV' in data.dtype.names:
                    uu_cov_data = data['UU_COV']
                else:
                    uu_cov_data = None
                if 'QQ_COV' in data.dtype.names:
                    qq_cov_data = data['QQ_COV']
                else:
                    qq_cov_data = None 
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
    For time-saving sake, the returned values are all in a single array of shape (2, 6, 50331648*).
    The 2 refers to the two frequencies of the fits files being used;
    The 6 refers to stokes_i, stokes_q, stokes_u, i_covariance, q_covariance, u_covariance, in that order;
    The 50331648 is the length of each of the arrays, this number is dependant on the specific fits file being input.

    Input: list/array/tuple of fits file names as strings, tuple of constants, array of frequencies of the files being input(externally created)
    Output: Array of shape (2, 6, 50331648*) containing stokes parameters and associated covariances at two different frequencies, all converted from units of Kcmb to MJy/sr
    """

    file_217, file_353 = files

    """ 
    The layout of the data array is:
    0,0: i_stokes_data_kcmb_217; 1,0: i_stokes_data_kcmb_353
    0,1: q_stokes_data_kcmb_217; 1,1: q_stokes_data_kcmb_353
    0,2: u_stokes_data_kcmb_217; 1,2: u_stokes_data_kcmb_353
    0,3: ii_cov_data_kcmb_217; 1,3: ii_cov_data_kcmb_353
    0,4: qq_cov_data_kcmb_217; 1,4: qq_cov_data_kcmb_353
    0,5: uu_cov_data_kcmb_217; 1,5: uu_cov_data_kcmb_353
    """

    extracted_data_array = np.array([
    extract_info(file_217)[:6],  # Take the first 6 elements of the output, ignore frequency and nest type
    extract_info(file_353)[:6]   # Take the first 6 elements of the output, ignore frequency and nest type
])

    conversion_factor_217, conversion_factor_353 = conversion_factors_Kcmb_to_MJy_sr(constants, frequencies)

    for i in range(3): #Stokes parameters are converted from Kcmb to MJy/sr
        extracted_data_array[0, i, :] = extracted_data_array[0, i, :] * conversion_factor_217
        extracted_data_array[1, i, :] = extracted_data_array[1, i, :] * conversion_factor_353
    for j in range(3,6): #Covariances are converted from Kcmb^2 to (MJy/sr)^2
        extracted_data_array[0, j, :] = extracted_data_array[0, j, :] * conversion_factor_217 ** 2
        extracted_data_array[1, j, :] = extracted_data_array[1, j, :] * conversion_factor_353 ** 2


    """ 
    The code past this comment decreases the resolution of the data from 2048 to 512 to decrease computing time. Delete past this comment and return the extracted_data_array if you want to keep the full resolution.
    """
    return extracted_data_array

def decrease_resolution(extracted_data_array, new_nside):
    extracted_data_array_new_nside = np.empty((2, 6, hp.nside2npix(new_nside)))

    # Fill the new array with downgraded maps
    for i in range(2):
        for j in range(6):
            extracted_data_array_new_nside[i, j, :] = hp.ud_grade(extracted_data_array[i, j, :], new_nside, order_in='NESTED', order_out='NESTED')


    return extracted_data_array_new_nside

def conversion_factors_Kcmb_to_MJy_sr(constants, frequencies):

    """ 
    This function computes the conversion factor between units of kcmb and MJy/sr as a function of frequency. The tuple of constants is called for this, and the last value, freq_0=353e9GHz is explicitely ignored.
    The equation used to convert between frequencies is the temperature derivative of the Planck function, it is a function of frequency.

    Inputs: Array of frequency values. These should be previously defined before calling this function. Tuple of constants.
    Returns: Array of conversion factors corresponding to the elements within the array of frequencies.
    """
    c, h, k, T_cmb, *ingore = constants

    T_cmb = 2.7255

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

def Chi2(parameters, frequencies, constants, stokes_arrays_correct_units, i):
    recreated_values = stokes_reconstruction(parameters, constants, frequencies)

    """
    This function quantifies the fit of the modelled emission vs the emission from the data. The recreated emission is calculated in the stokes_reconstruction function.
    Chi^2 values are calculated for each of the 3 stokes parameters at each of the 2 input frequencies, then the sum of all of these 6 values is taken and returned as the Chi^2 value.
    This recreated emission values are a function of the 'parameters' array that is called, the indices of this array will be optimized to minimize the value returned by this function.
    This function calculates the Chi^2 value for a single pixel, it must be used in a loop over all pixels of the data that you are trying to recreate.
    The stokes_arrays_correct_units element that is called is the data from the fits file after being converted from Kcmb to MJy/sr. It must be extracted and called from outside of the function to greatly improve efiiciency.

    Inputs: Parameters array(will be optimized), frequencies, constants, stokes_arrays_correct_units, i (i is the index of the loop that this function must be run through)
    Outputs: Sum of all Chi^2 values for each stokes parameter across the array of frequencies
    """

    #To optimize all three parameters, use this set:
    Chi2_Stokes_I = (recreated_values[:, 0] - stokes_arrays_correct_units[:, 0, i]) ** 2 / stokes_arrays_correct_units[:, 3, i]
    Chi2_Stokes_Q = (recreated_values[:, 1] - stokes_arrays_correct_units[:, 1, i]) ** 2 / stokes_arrays_correct_units[:, 4, i]
    Chi2_Stokes_U = (recreated_values[:, 2] - stokes_arrays_correct_units[:, 2, i]) ** 2 / stokes_arrays_correct_units[:, 5, i]

    #To ignore the optimization of a parameter, uncomment the following lines and comment out the above lines:
    #Chi2_Stokes_I = np.zeros(2)
    #Chi2_Stokes_Q = np.zeros(2)
    #Chi2_Stokes_U = np.zeros(2)

    return np.sum(Chi2_Stokes_I) + np.sum(Chi2_Stokes_Q) + np.sum(Chi2_Stokes_U)

def define_bounds():
    return [(-np.infty, np.infty), #T
          (-np.infty, np.infty), #Beta
          (-np.infty, np.infty), #Tau
          (-np.infty, np.infty), #Psi(radians)
          (-np.infty, np.infty), #Alpha
          (-np.infty, np.infty)] #p_frac