# create modular code for creating vEPAC pulse from measurement or numerically
import yaml
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, c
from lasy.optical_elements import PolynomialSpectralPhase

from lasy.profiles.combined_profile import CombinedLongitudinalTransverseProfile
from lasy.profiles import FromInsightFile
from lasy.laser import Laser

import profiles_from_measurement as pfm
import profiles_from_parameters as pfp


# open config file to extract needed parameters
with open("C:\\Users\\skq54543\\OneDrive - Science and Technology Facilities Council\\EPAC Simulation\\actual_modules_and_script\\config.yaml", "r") as f:
    config = yaml.full_load(f)

laser_params = config['laser_params']
measurement_params = config['measurement_params']
analytic_params = config['analytic_params']
grid_params = config['grid_params']
figure_params = config['figure_params']

energy_J = laser_params['energy_J']
wavelength = laser_params['wavelength']
polarisation = laser_params['polarisation']

dropout_figs = figure_params['dropout_figs']

no_spatial_steps = grid_params['no_spatial_steps'] # in radial direction
no_time_steps = grid_params['no_time_steps']
n_azimuthal_modes = grid_params['n_azimuthal_modes']
n_theta_evals = grid_params['n_theta_evals'] # default is 2*n_azimuthal_modes - 1

near_field_image = measurement_params['is_NF_measurement']
transverse_filepath = measurement_params['NF_filepath']
cal = measurement_params['pixel_size']
spectrum_measurement = measurement_params['is_FROG_measurement']
longitudinal_filepath = measurement_params['FROG_filepath']
insight = measurement_params['is_INSIGHT']
insight_filepath = measurement_params['INSIGHT_filepath']

beam_waist = analytic_params['beam_waist']
f_number = pi * beam_waist / (2*wavelength)
na = 1 / (2*f_number)
pixel_scale = 10 * beam_waist / no_spatial_steps
coherence_factor = analytic_params['coherence_factor']
is_zernike = analytic_params['zernike']
zernike_dict = {tuple(map(int, k.split(','))): v for k, v in analytic_params["zernike_dict"].items()}

pulse_duration = float(analytic_params['pulse_duration'])
t_peak = analytic_params['t_peak']
cep_phase = analytic_params['cep_phase']
spectral_phase = np.zeros(no_time_steps) # spectral phase profile in frequency domain (spectral phase of a Gaussian is zero)
gdd = analytic_params['gdd']
tod = analytic_params['tod']
fod = analytic_params['fod']


if insight:
    laser_profile = pfm.ProfileFromInsight(insight_filepath, polarisation, wavelength, dropout_figs=dropout_figs)
else:
    # for now defining numeric pulse as Gaussian and adding aberrations and spectral phase - could change to supergaussian etc later on
    if near_field_image:
        # define transverse spectrum from image (import module)
        transverse_profile, measured_waist = pfm.TransverseProfileFromMeasurement(transverse_filepath, cal, dropout_figs=dropout_figs)
    else:
        spatial_steps = np.linspace(-5*beam_waist, 5*beam_waist, no_spatial_steps) # could define this the same every time instead of from inputs
        x_coords, y_coords = np.meshgrid(spatial_steps, spatial_steps)
        transverse_profile = pfp.TransverseProfileFromParameters(x_coords, y_coords, beam_waist, wavelength, pixel_scale, na, no_spatial_steps, coherence_factor, dropout_figs=dropout_figs, is_zernike=is_zernike, zernike_dict=zernike_dict)

    if spectrum_measurement:
        # define longitudinal spectrum from FROG measurement
        longitudinal_profile = pfm.LongitudinalProfileFromMeasurement(longitudinal_filepath, dropout_figs=dropout_figs)
    else:
        # add functionality for time step
        time_steps = np.linspace(-5*pulse_duration, 5*pulse_duration, no_time_steps) # could define this the same every time instead of from inputs
        longitudinal_profile = pfp.LongitudinalProfileFromParameters(time_steps, t_peak, pulse_duration, wavelength, beam_waist, cep_phase, gdd=gdd, tod=tod, fod=fod, dropout_figs=dropout_figs)

    # combine longitudinal and transverse profiles into single profile
    laser_profile = CombinedLongitudinalTransverseProfile(
        wavelength=wavelength,
        pol=polarisation,
        laser_energy=energy_J,
        long_profile=longitudinal_profile,
        trans_profile=transverse_profile,
    )



# Create a LASY Laser object
# Cylindrical geometry easiest for interfacing with FBPIC

dimensions = "rt"  # Use cylindrical geometry (better for fbpic)
num_points = (no_spatial_steps, no_time_steps) 


# choosing lo and hi based on pulse
if near_field_image:
    radial_limit = 5*measured_waist
else:
    radial_limit = 30*beam_waist


if spectrum_measurement:
    time_limit = 100e-15 # decide on 100fs domain for now (same domain as used in LongitudinalProfileFromMeasurement) maybe change this later depending on what pulses we think will get
else:
    time_limit = 5*pulse_duration


lo = (0, -time_limit)
hi = (radial_limit, time_limit)


laser = Laser(dimensions, lo, hi, num_points, laser_profile, n_azimuthal_modes=n_azimuthal_modes)

# add gdd, tod, fod if defining profile from parameters
if spectrum_measurement == False:
    omega0 = 2*pi*c / wavelength
    dazzler = PolynomialSpectralPhase(omega0, gdd, tod, fod)
    laser.apply_optics(dazzler)

# Laser object as drop out figure?
laser.show()
plt.title("Laser Object")
plt.show()

# write laser object to file
file_prefix = "test_pulse_from_parameters"  # The file name will start with this prefix
file_format = "h5"  # Format to be used for the output file
directory = 'C:\\Users\\skq54543\\OneDrive - Science and Technology Facilities Council\\EPAC Simulation' # CHANGE THIS!

laser.write_to_file(file_prefix, file_format, write_dir='C:\\Users\\skq54543\\OneDrive - Science and Technology Facilities Council\\EPAC Simulation')

