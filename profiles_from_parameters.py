# initialise transverse and longitudinal profiles from given parameters
import astropy.units as u
import pyoptica as po
import numpy as np
from scipy.constants import pi, c
import matplotlib.pyplot as plt

from lasy.profiles.longitudinal import LongitudinalProfileFromData
from lasy.profiles.transverse import TransverseProfileFromData

def transverse_gaussian(x_coords, y_coords, beam_waist):
    return(np.exp(-(x_coords**2 + y_coords**2) / beam_waist**2))

def longitudinal_gaussian(t, t_peak, pulse_duration, beam_waist, cep_phase):
    return np.exp(-((t - t_peak)**2) / pulse_duration**2 + 1.0j *(cep_phase + beam_waist*t_peak))

# need to initialise grid to calculate profiles on - do in this module or in main script?


def TransverseProfileFromParameters(x_coords, y_coords, beam_waist, wavelength, pixel_scale, na, npix, coherence_factor, dropout_figs=True, is_zernike=False, zernike_dict=None):
    # define transverse gaussian from grid and beam waist
    # propagate using PyOptica to add Zernike aberrations
    # inputs: beam waist (f_number), dictionary of Zernike aberrations
    
    transverse_spectrum = transverse_gaussian(x_coords, y_coords, beam_waist)

    # IMAGE WITH ARBITRARY ZERNIKE POLYNOMIAL ABERRATIONS IN THE BEAM
    # NEED TO REMEMBER TO CONVERT TO ASTROPY UNITS WHEN USING PYOPTICA
    wavelength_unit = wavelength * u.m
    pixel_scale_unit = pixel_scale * u.m

    img_system = po.ImagingSystem(wavelength_unit, pixel_scale_unit, npix, na=na, coherence_factor=coherence_factor)

    if is_zernike:
        img_system.load_zernikes(zernike_dict, 'mn') # remember about the convention

    img_system.calculate()  # We need to run this to precalculate the pupil!

    wf = po.Wavefront(wavelength_unit, pixel_scale_unit, npix)
    wf.intensity = transverse_spectrum

    # imaging with aberrations 
    image = img_system.image_wavefront(wf)
    image_array = image.image

    transverse_profile = TransverseProfileFromData(image_array, lo=[-5*beam_waist,-5*beam_waist], hi=[5*beam_waist,5*beam_waist])

    if dropout_figs:

        _ = image.plot(image=dict(vmax=0.3), fig_options=dict(dpi=130))
        plt.title('Transverse profile (from parameters)')
        plt.xlabel ('x (m)')
        plt.ylabel('y (m)')

        fig_caption = "Beam waist = " + str(beam_waist) + "m\nZernikes {(m, n): value} = " + str(zernike_dict) + "\n(Gaussian)"
        plt.figtext(0.65, 0, fig_caption, horizontalalignment='center')

        plt.gcf().subplots_adjust(bottom=0.15)
        plt.show()

    return transverse_profile


def spectral_phase_calculator_frequency(omega, central_wavelength, total_phase=0, group_delay=0, gdd=0, tod=0, fod=0):
    omega0 = 2 * pi * c / central_wavelength
    omega_difference = omega - omega0
    spectral_phase = total_phase + group_delay*(omega_difference) +1/2*gdd*(omega_difference)**2 + 1/6*tod*(omega_difference)**3 + 1/24*fod*(omega_difference)**4
    return spectral_phase

def LongitudinalProfileFromParameters(time_steps, t_peak, pulse_duration, wavelength, beam_waist, cep_phase, gdd=0, tod=0, fod=0, dropout_figs=True):
    # define longitudinal gaussian from grid and pulse duration
    # annoyingly have to initialise laser pulse before you can add gdd, tod, etc (do this in main script)
   
    longitudinal_spectrum = longitudinal_gaussian(time_steps, t_peak, pulse_duration, beam_waist, cep_phase)
    spectral_phase = np.zeros(len(time_steps))
    temporal_phase = np.fft.ifft(spectral_phase)

    longitudinal_data = {
    "datatype": "temporal",
    "wavelength": wavelength,
    "axis": time_steps, # time axis of pulse duration measurement
    "intensity": longitudinal_spectrum, # vertical axis of pulse duration measurement
    "dt": 2*pulse_duration / len(time_steps),
    "phase": temporal_phase,
    }

    longitudinal_profile = LongitudinalProfileFromData(longitudinal_data, lo=-5*pulse_duration, hi=5*pulse_duration)

    if dropout_figs:

        omega_min = 2*pi*c / (wavelength+100e-9)
        omega_max =  2*pi*c / (wavelength-100e-9)
        no_points = 500

        omega_array = np.linspace(omega_min, omega_max, no_points)

        spectral_phase_profile = spectral_phase_calculator_frequency(omega_array, wavelength, total_phase=0, gdd=gdd, tod=tod, fod=fod)

        fig, ax = plt.subplots(1, 2, figsize=(12, 4), tight_layout=True)

        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Intensity") # UNIT? (a.u.)
        ax[0].set_title('Temporal intensity (Longitudinal Profile)')
        ax[0].plot(time_steps, longitudinal_spectrum, color='blue')

        fig_caption_1 = "Pulse duration = " + str(pulse_duration) + "s\n(Gaussian)"
        ax[0].annotate(fig_caption_1, xy=(0.25, 0.05), xytext = (0.25, 0.05), xycoords='figure fraction', horizontalalignment='center')


        ax[1].set_xlabel("Frequency (Hz)")
        ax[1].set_ylabel("Spectral Phase") #UNIT? AU
        ax[1].set_title('Spectral phase of laser pulse (Gaussian + gdd,tod,fod)')
        ax[1].plot(omega_array, spectral_phase_profile, color='red')

        fig_caption_2 = "GDD = " + str(gdd) + "fs^2\nTOD = " + str(tod) + "fs^3\nFOD = " + str(fod) + "fs^4"
        ax[1].annotate(fig_caption_2, xy=(0.75, 0), xytext = (0.75,0), xycoords='figure fraction', horizontalalignment='center')

        plt.gcf().subplots_adjust(bottom=0.2)
        plt.tight_layout()

        plt.show()

    return longitudinal_profile