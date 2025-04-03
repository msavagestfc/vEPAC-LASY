import yaml

# all values should be given in SI units

laser_params = {
    'energy_J': 1,
    'wavelength': 800e-9,
    'polarisation': (1,0),
}

measurement_params = {
    'is_NF_measurement':False,
    'NF_filepath':"C:\\Users\\skq54543\\OneDrive - Science and Technology Facilities Council\\EPAC Simulation\\transverse_profile_test.png",
    'pixel_size':0.2e-6,
    'is_FROG_measurement':True,
    'FROG_filepath':"https://github.com/user-attachments/files/17414077/df_intensity_spectral_v3.csv",
    'is_INSIGHT': False,
    'INSIGHT_filepath': "C:\\Users\\skq54543\\OneDrive - Science and Technology Facilities Council\\EPAC Simulation\\INSIGHT measurements\\Exyt_0_test.h5",
}

analytic_params = {
    'beam_waist': 5e-6,
    'pulse_duration': 30e-15,
    't_peak':0.0,  # Location of the peak of the laser pulse in time
    'cep_phase':0, # carrier envelope phase
    'n_azimuthal_modes': 3,
    'n_theta_evals': 5,
    'num_points': (500, 500), # (r,t)
    'coherence_factor': 0,
    'zernike_dict': {(1,3):0.1, (2,2):0.1},
    'gdd': 21000e-30,
    'tod': -20000e-45,
    'fod': 800000e-60,
}

grid_params = {
    'npix': 500, # no.of pixels for imaging through Zernike optics - used to initialise spatial frid
    'no_time_steps': 500, # used to initialise time grid
}

figure_params = {
    'dropout_figs':True,
}

data = {"laser_params":laser_params, "measurement_params":measurement_params, "analytic_params":analytic_params,
        "grid_params": grid_params, 'figure_params': figure_params}

# write to yaml file config.yaml
with open("config.yaml", "w") as f:
    for key, value in data.items():
        yaml.dump({key: value}, f, default_flow_style=False)
        f.write("\n\n")