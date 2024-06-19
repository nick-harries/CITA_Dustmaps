import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from statistics import mode


def plot_recreated_values_and_fit(stokes_arrays_reconstructed, stokes_arrays_correct_units, arrays_to_optimize, starting_index, ending_index, frequencies):
    """
    This function calls the array created by the execute_Chi2_optimization function, as well as the Planck data. It then plots the modelled emission next to the true emission for each frequency and stokes parameter.
    """
    Stokes_parameters = ['I', 'Q', 'U']
    for freq in range(np.shape(stokes_arrays_reconstructed)[0]):
        for parameter in range(np.shape(stokes_arrays_reconstructed)[1]):
            plt.figure()
            if [freq, parameter] in arrays_to_optimize:                
                residual = stokes_arrays_correct_units[freq, parameter, starting_index:ending_index] - stokes_arrays_reconstructed[freq, parameter, :]
                plt.subplot(2,2 ,1)
                plt.title('Stokes %s %.2e (boundless)'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot(stokes_arrays_correct_units[freq, parameter, starting_index:ending_index], color='red', label='data')
                plt.plot(stokes_arrays_reconstructed[freq, parameter], color='blue', label='model')
                plt.xlabel('Index Number')
                plt.ylabel('Intensity (MJy/sr)')
                plt.legend(loc='upper right')
                
                plt.subplot(2, 2, 2)                
                plt.title('Stokes %s %.2e data'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot(stokes_arrays_correct_units[freq, parameter, starting_index:ending_index], color='red', label='data')
                plt.xlabel('Index Number')
                plt.ylabel('Intensity (MJy/sr)')

                plt.subplot(2, 2, 4)
                plt.title('Stokes %s %.2e residual'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot(residual, color='black')
                plt.xlabel('Index Number')
                plt.ylabel('Intensity (MJy/sr)')
            else:
                plt.subplot(2, 2, 1)
                plt.title('No data Stokes %s %.2e'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot()

                plt.subplot(2, 2, 2)
                plt.title('No data Stokes %s %.2e'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot()

                plt.subplot(2, 2, 4)
                plt.title('No data Stokes %s %.2e'%(Stokes_parameters[parameter], frequencies[freq]))
                plt.plot()

            plt.subplot(2, 2, 3)
            plt.title('Stokes %s %.2e GHz model'%(Stokes_parameters[parameter], frequencies[freq]))
            plt.plot(stokes_arrays_reconstructed[freq, parameter], color='blue')
            plt.xlabel('Index Number')
            plt.ylabel('Intensity (MJy/sr')
            plt.tight_layout()
        plt.show()
                
def plot_optimized_parameters(optimized_parameters_array):
    """
    This function calls the array created by the  execute_Chi2_optimization function and the plots the optimized parameters as a function of index. These plots are titled with which parameter they show as well as the mean value across the sample.
    """
    Temperature, Beta, Tau, Psi, Alpha, p_frac= data[:, -7], data[:, -6], data[:, -5], data[:, -4], data[:, -3], data[:, -2]

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
    Temperature, Beta, Tau, Psi, Alpha, p_frac= data[:, -7], data[:, -6], data[:, -5], data[:, -4], data[:, -3], data[:, -2]

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

def plot_minimized_Chi2(optimized_parameters_array):
    """
    This function calls the array created by the execute_Chi2_optimization function. It then plots the Chi^2 value as a function of pixel index and prints the mean Chi^2 across the sample in the title.
    """
    plt.figure()
    plt.title('Chi^2 value by pixel, mean: %.2f'%np.mean(optimized_parameters_array[-1]))
    plt.plot(optimized_parameters_array[-1])
    plt.xlabel('Index number')
    plt.ylabel('Chi^2 value')
   #plt.savefig('Minimized_Chi2.png')
    plt.show()
