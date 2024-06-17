import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from statistics import mode


def plot_recreated_values_and_fit(data, stokes_arrays_correct_units, starting_index, ending_index):
    """
    This function calls the array created by the execute_Chi2_optimization function, as well as the Planck data. It then plots the modelled emission next to the true emission for each frequency and stokes parameter.
    """
    stokes_i_217_recreated, stokes_i_353_recreated = data[:, 0], data[:, 3]
    stokes_q_217_recreated, stokes_q_353_recreated = data[:, 1], data[:, 4]
    stokes_u_217_recreated, stokes_u_353_recreated = data[:, 2], data[:, 5]

    stokes_i_217_data, stokes_i_353_data = stokes_arrays_correct_units[0, 0, starting_index:ending_index], stokes_arrays_correct_units[1, 0, starting_index:ending_index]
    stokes_q_217_data, stokes_q_353_data = stokes_arrays_correct_units[0, 1, starting_index:ending_index], stokes_arrays_correct_units[1, 1, starting_index:ending_index]
    stokes_u_217_data, stokes_u_353_data = stokes_arrays_correct_units[0, 2, starting_index:ending_index], stokes_arrays_correct_units[1, 2, starting_index:ending_index]

    i_217_residual = stokes_i_217_data - stokes_i_217_recreated
    i_353_residual = stokes_i_353_data - stokes_i_353_recreated

    q_217_residual = stokes_q_217_data - stokes_q_217_recreated
    q_353_residual = stokes_q_353_data - stokes_q_353_recreated

    u_217_residual = stokes_u_217_data - stokes_u_217_recreated
    u_353_residual = stokes_u_353_data - stokes_u_353_recreated


    plt.figure()
    plt.subplot(221)
    plt.title('Stokes I 217 (boundless)')
    plt.plot(stokes_i_217_data, color='red', label = 'data')
    plt.plot(stokes_i_217_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data I 217GHz')
    plt.plot(stokes_i_217_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('I model (boundless)')
    plt.plot(stokes_i_217_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('I 217GHz residual (data-model)')
    plt.plot(i_217_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(221)
    plt.title('Stokes I 353 (boundless)')
    plt.plot(stokes_i_353_data, color='red', label = 'data')
    plt.plot(stokes_i_353_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data I 353GHz')
    plt.plot(stokes_i_353_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('I model (boundless)')
    plt.plot(stokes_i_353_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('I 353GHz residual (data-model)')
    plt.plot(i_353_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(221)
    plt.title('Stokes Q 217 (boundless)')
    plt.plot(stokes_q_217_data, color='red', label = 'data')
    plt.plot(stokes_q_217_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data Q 217GHz')
    plt.plot(stokes_q_217_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('Q model (boundless)')
    plt.plot(stokes_q_217_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('Q 217GHz residual (data-model)')
    plt.plot(q_217_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(221)
    plt.title('Stokes Q 353 (boundless)')
    plt.plot(stokes_q_353_data, color='red', label = 'data')
    plt.plot(stokes_q_353_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data Q 353GHz')
    plt.plot(stokes_q_353_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('Q model (boundless)')
    plt.plot(stokes_q_353_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('Q 353GHz residual (data-model)')
    plt.plot(q_353_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(221)
    plt.title('Stokes U 217 (boundless)')
    plt.plot(stokes_u_217_data, color='red', label = 'data')
    plt.plot(stokes_u_217_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data U 217GHz')
    plt.plot(stokes_u_217_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('U model (boundless)')
    plt.plot(stokes_u_217_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('U 217GHz residual (data-model)')
    plt.plot(u_217_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.subplot(221)
    plt.title('Stokes U 353 (boundless)')
    plt.plot(stokes_u_353_data, color='red', label = 'data')
    plt.plot(stokes_u_353_recreated, color='blue', label = 'recreation')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(222)
    plt.title('Raw Planck data U 353GHz')
    plt.plot(stokes_u_353_data, color='red', label = 'data')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(223)
    plt.title('U model (boundless)')
    plt.plot(stokes_u_353_recreated, color='blue', label='model')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')

    plt.subplot(224)
    plt.title('U 353GHz residual (data-model)')
    plt.plot(u_353_residual, label='residual')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

def plot_optimized_parameters(data):
    """
    This function calls the array created by the  execute_Chi2_optimization function and the plots the optimized parameters as a function of index. These plots are titled with which parameter they show as well as the mean value across the sample.
    """
    Temperature, Beta, Tau, Psi, Alpha, p_frac= data[:, 6], data[:, 7], data[:, 8], data[:, 9], data[:, 10], data[:, 11]

    #plt.subplot(231)
    plt.plot(Temperature)
    plt.title('Temperature, mean: %.2f, mode: %.2f'%(np.mean(Temperature), mode(Temperature)))
    plt.xlabel('Index Number')
    plt.ylabel('Magnitude in Kelvin')
   #plt.savefig('Temperature_array.png')

    #plt.subplot(232)
    plt.figure()
    plt.plot(Beta)
    plt.title('beta, mean: %.2f, mode: %.2f'%(np.mean(Beta), mode(Temperature)))
    plt.xlabel('Index Number')
    plt.ylabel('Unitless Magnitude')
   #plt.savefig('Beta_array.png')

    #plt.subplot(233)
    plt.figure()
    plt.plot(Tau)
    plt.title('tau, mean: %.2e, mode: %.2f'%(np.mean(Tau), mode(Tau)))
    plt.xlabel('Index Number')
    plt.ylabel('Unitless Magnitude')
   #plt.savefig('Tau_array.png')

    #plt.subplot(234)
    plt.figure()
    plt.plot(Psi)
    plt.title('psi, mean: %.2f, mode: %.2f'%(np.mean(Psi), mode(Psi)))
    plt.xlabel('Index Number')
    plt.ylabel('Magnitude in radians')
   #plt.savefig('Psi_array.png')

    #plt.subplot(235)
    plt.figure()
    plt.plot(Alpha)
    plt.title('alpha, mean: %.2f, mode: %.2f'%(np.mean(Alpha), mode(Alpha)))
    plt.xlabel('Index Number')
    plt.ylabel('Magnitude in radians')
   #plt.savefig('Alpha_array.png')

    #plt.subplot(236)
    plt.figure()
    plt.plot(p_frac)
    plt.title('p_frac, mean: %.2f, mode: %.2f'%(np.mean(p_frac), mode(p_frac)))
    plt.xlabel('Index Number')
    plt.ylabel('Fractional Polarization')
   #plt.savefig('Fractional_polarization_array.png')

    #plt.tight_layout()

    plt.show()

def plot_optimized_parameters_histograms(data):
    """
    This function calls the array created by the  execute_Chi2_optimization function and the plots the optimized parameters as a function of index. These plots are titled with which parameter they show as well as the mean value across the sample.
    """
    Temperature, Beta, Tau, Psi, Alpha, p_frac = data[:, 6], data[:, 7], data[:, 8], data[:, 9], data[:, 10], data[:, 11]

    #plt.subplot(231)
    plt.figure()
    plt.hist(Temperature, bins = 100)
    plt.title('Temperature, mean: %.2f, mode : %.2f'%(np.mean(Temperature), mode(Temperature)))
    plt.xlabel('Degrees Kelvin')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Temperature_histogram_Q_only_boundless.png')

    #plt.subplot(232)
    plt.figure()
    plt.hist(Beta, bins = 100)
    plt.title('beta, mean: %.2f, mode: %.2f'%(np.mean(Beta), mode(Beta)))
    plt.xlabel('Beta Value')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Beta_histogram_Q_only_boundless.png')

    #plt.subplot(233)
    plt.figure()
    plt.hist(Tau, bins = 100)
    plt.title('tau, mean: %.2e, mode: %.2e'%(np.mean(Tau), mode(Tau)))
    plt.xlabel('Tau Value')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Tau_histogram_Q_only_boundless.png')

    #plt.subplot(234)
    plt.figure()
    plt.hist(Psi, bins = 100)
    plt.title('psi, mean: %.2f, mode: %.2f'%(np.mean(Psi), mode(Psi)))
    plt.xlabel('Psi value (radians)')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Psi_histogram_Q_only_boundless.png')

    #plt.subplot(235)
    plt.figure()
    plt.hist(Alpha, bins = 100)
    plt.title('alpha, mean: %.2f, mode: %.2f'%(np.mean(Alpha), mode(Alpha)))
    plt.xlabel('Alpha value (radians)')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Alpha_histogram_Q_only_boundless.png')

    #plt.subplot(236)
    plt.figure()
    plt.hist(p_frac, bins = 100)
    plt.title('p_frac, mean: %.2f, mode: %.2f'%(np.mean(p_frac), mode(p_frac)))
    plt.xlabel('Fractional Polarization values')
    plt.ylabel('Number of Occurrences')
    #plt.savefig('Fractional_polarization_histogram_Q_only_boundless.png')

    #plt.tight_layout()
    #plt.savefig('Optimized_parameters_histograms.png')
    plt.show()

    plt.show()

def plot_raw_data_against_itself(stokes_arrays_correct_units, starting_index, ending_index):
    """
    This function plots the data within the range that has been sampled. It plots the I, Q, and U from a common frequency on the same plot. It does this for both stokes parameters value as well as covariance value. 
    The plots are made with logarithmic Y-axes to show differences of intensities between parameters.
    """
    stokes_i_217, stokes_q_217, stokes_u_217 = stokes_arrays_correct_units[0, 0, starting_index:ending_index], stokes_arrays_correct_units[0, 1, starting_index:ending_index], stokes_arrays_correct_units[0, 2, starting_index:ending_index]
    stokes_i_353, stokes_q_353, stokes_u_353 = stokes_arrays_correct_units[1, 0, starting_index:ending_index], stokes_arrays_correct_units[1, 1, starting_index:ending_index], stokes_arrays_correct_units[1, 2, starting_index:ending_index]

    ii_cov_217, qq_cov_217, uu_cov_217 = stokes_arrays_correct_units[0, 3, starting_index:ending_index], stokes_arrays_correct_units[0, 4, starting_index:ending_index], stokes_arrays_correct_units[0, 5, starting_index:ending_index]
    ii_cov_353, qq_cov_353, uu_cov_353 = stokes_arrays_correct_units[1, 3, starting_index:ending_index], stokes_arrays_correct_units[1, 4, starting_index:ending_index], stokes_arrays_correct_units[1, 5, starting_index:ending_index]

    plt.figure()
    plt.title('Stokes 217')
    plt.plot(stokes_i_217, label='I')
    plt.plot(stokes_q_217, label='Q')
    plt.plot(stokes_u_217, label='U')
    plt.legend(loc='upper right')
    plt.yscale('log')

    plt.figure()
    plt.title('Stokes 353')
    plt.plot(stokes_i_353, label='I')
    plt.plot(stokes_q_353, label='Q')
    plt.plot(stokes_u_353, label='U')
    plt.legend(loc='upper right')
    plt.yscale('log')

    plt.figure()
    plt.title('Covariances 217')
    plt.plot(ii_cov_217, label='I')
    plt.plot(qq_cov_217, label='Q')
    plt.plot(uu_cov_217, label='U')
    plt.legend(loc='upper right')
    plt.yscale('log')

    plt.figure()
    plt.title('Covariances 353')
    plt.plot(ii_cov_353, label='I')
    plt.plot(qq_cov_353, label='Q')
    plt.plot(uu_cov_353, label='U')
    plt.legend(loc='upper right')
    plt.yscale('log')

    plt.show()

def plot_raw_data(stokes_arrays_correct_units, starting_index, ending_index):
    """
    This function plots the data within the range that has been sampled. It plots the I, Q, and U from a common frequency on the same plot. It does this for both stokes parameters value as well as covariance value. 
    The plots are made with logarithmic Y-axes to show differences of intensities between parameters.
    """
    stokes_i_217, stokes_q_217, stokes_u_217 = stokes_arrays_correct_units[0, 0, starting_index:ending_index], stokes_arrays_correct_units[0, 1, starting_index:ending_index], stokes_arrays_correct_units[0, 2, starting_index:ending_index]
    stokes_i_353, stokes_q_353, stokes_u_353 = stokes_arrays_correct_units[1, 0, starting_index:ending_index], stokes_arrays_correct_units[1, 1, starting_index:ending_index], stokes_arrays_correct_units[1, 2, starting_index:ending_index]

    ii_cov_217, qq_cov_217, uu_cov_217 = stokes_arrays_correct_units[0, 3, starting_index:ending_index], stokes_arrays_correct_units[0, 4, starting_index:ending_index], stokes_arrays_correct_units[0, 5, starting_index:ending_index]
    ii_cov_353, qq_cov_353, uu_cov_353 = stokes_arrays_correct_units[1, 3, starting_index:ending_index], stokes_arrays_correct_units[1, 4, starting_index:ending_index], stokes_arrays_correct_units[1, 5, starting_index:ending_index]

    plt.figure()
    plt.subplot(211)
    plt.title('Stokes I 217')
    plt.plot(stokes_i_217, label='I 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('Stokes I 353')
    plt.plot(stokes_i_353, label='I 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('stokes_i_raw.png')

    plt.figure()
    plt.subplot(211)
    plt.title('Stokes U 217')
    plt.plot(stokes_u_217, label='U 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('Stokes U 353')
    plt.plot(stokes_u_353, label='U 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('stokes_u_raw.png')

    plt.figure()
    plt.subplot(211)
    plt.title('Stokes Q 217')
    plt.plot(stokes_q_217, label='q 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('Stokes Q 353')
    plt.plot(stokes_q_353, label='Q 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('stokes_q_raw.png')

    plt.figure()
    plt.subplot(211)
    plt.title('I Covariance 217')
    plt.plot(ii_cov_217, label='I Covariance 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('I Covariance 353')
    plt.plot(ii_cov_353, label='I Covariance 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)^2')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('i_cov_raw.png')

    plt.figure()
    plt.subplot(211)
    plt.title('Q Covariance 217')
    plt.plot(qq_cov_217, label='Q Covariance 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('Q Covariance 353')
    plt.plot(qq_cov_353, label='Q Covariance 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)^2')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('q_cov_raw.png')


    plt.figure()
    plt.subplot(211)
    plt.title('U Covariance 217')
    plt.plot(uu_cov_217, label='U Covariance 217GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)')
    plt.legend(loc='upper right')
    plt.subplot(212)
    plt.title('U Covariance 353')
    plt.plot(uu_cov_353, label='U Covariance 353GHz')
    plt.xlabel('Index Number')
    plt.ylabel('Intensity (MJy/sr)^2')
    plt.legend(loc='upper right')
    plt.tight_layout()
   #plt.savefig('u_cov_raw.png')

    plt.show()

def plot_minimized_Chi2(data):
    """
    This function calls the array created by the execute_Chi2_optimization function. It then plots the Chi^2 value as a function of pixel index and prints the mean Chi^2 across the sample in the title.
    """
    plt.figure()
    plt.title('Chi^2 value by pixel, mean: %.2f'%np.mean(data[:, 12]))
    plt.plot(data[:, 12])
    plt.xlabel('Index number')
    plt.ylabel('Chi^2 value')
   #plt.savefig('Minimized_Chi2.png')
    plt.show()