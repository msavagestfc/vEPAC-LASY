import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from PIL import Image

from lasy.profiles.longitudinal import LongitudinalProfileFromData
from lasy.profiles.transverse import TransverseProfileFromData
from lasy.profiles import FromInsightFile

from lasy.utils.mode_decomposition import estimate_best_HG_waist

def LongitudinalProfileFromMeasurement(longitudinal_filepath, dropout_figs=True, lo=-100e-15, hi=100e-15):
# take in filepath from FROG measurement and return longitudinal laser profile

    exp_frequency = np.loadtxt(longitudinal_filepath, usecols=0, dtype="float")  # Hz
    exp_spectrum = np.loadtxt(longitudinal_filepath, usecols=1, dtype="float")  # Arbitary units
    exp_phase = np.loadtxt(longitudinal_filepath, usecols=2, dtype="float")  # rad

    # Initialise a LASY LongitudinalProfile from experimentally measured spectrum
    # central wavelength is calculated in this step and used later on
    longitudinal_data = {
        "datatype": "spectral",
        "axis_is_wavelength": False,
        "axis": exp_frequency,
        "intensity": exp_spectrum,
        "phase": exp_phase,
        "dt": 1e-15, # femtosecond time step (can change this later)
    }

    # Create the longitudinal profile. The temporal range is from -200 to +200 femtoseconds
    # HAD TO CHANGE LASY SOURCE CODE TO GET THIS FUNCTION TO WORK
    longitudinal_profile = LongitudinalProfileFromData(
        longitudinal_data, lo=lo, hi=hi
    )

    if dropout_figs:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4), tight_layout=True)

        # Spectral data
        exp_spectrum /= np.max(exp_spectrum)  # Normalize the spectrum
        color = "tab:red"
        axs[0].set_xlabel("Frequency (Hz)", fontsize=12)
        axs[0].set_ylabel("Spectral intensity (AU)", color=color, fontsize=12)
        axs[0].plot(exp_frequency, exp_spectrum, color=color)
        axs[0].tick_params(axis="y", labelcolor=color)
        axs[0].set_title("Measured data (Spectral space)", fontsize=15)
        ax0 = axs[0].twinx()
        color = "tab:blue"
        ax0.set_ylabel("Spectral phase (radian)", color=color, fontsize=12)
        ax0.plot(exp_frequency, exp_phase, color=color)
        ax0.tick_params(axis="y", labelcolor=color)


        # Temporal data
        color = "tab:red"
        axs[1].set_xlabel("Time (s)", fontsize=12)
        axs[1].set_ylabel("Amplitude (AU)", color=color, fontsize=12)
        axs[1].plot(
            longitudinal_profile.time,
            np.sqrt(longitudinal_profile.temporal_intensity),
            color=color,
        )
        axs[1].tick_params(axis="y", labelcolor=color)
        axs[1].set_title("Reconstructed data (Temporal space)", fontsize=15)
        ax1 = axs[1].twinx()
        color = "tab:blue"
        ax1.set_ylabel("Temporal phase (radian)", color=color, fontsize=12)
        ax1.plot(longitudinal_profile.time, longitudinal_profile.temporal_phase, color=color)
        ax1.tick_params(axis="y", labelcolor=color)

        plt.show()


    # return profile and upper limit of time domain (for initialising laser object)
    return longitudinal_profile


def TransverseProfileFromMeasurement(transverse_filepath, pixel_size, dropout_figs=True):
# take in near field image and return transverse laser profile
    img = Image.open(transverse_filepath)
    intensity_data = np.array(img)
    intensity_scale = np.max(intensity_data)  # Maximum value of the intensity

    # ADD FUNCTIONALITY FOR CLEANING IMAGE
    intensity_data[intensity_data < intensity_scale / 100] = 0 # very simple cleaning method
    # or could use Hermite-Gauss decomposition

    # create a LASY TransverseProfile from experimental data
    nx, ny = intensity_data.shape
    lo = (0, 0)  # Lower bounds in x and y
    hi = (ny * pixel_size, nx * pixel_size)  # Upper bounds in x and y

    # Create the transverse profile. This also centers the data by default
    # THIS IS NOT CENTERING THE DATA BY DEFAULT!!!
    transverse_profile = TransverseProfileFromData(
        intensity_data, [lo[0], lo[1]], [hi[0], hi[1]], center_data = True,
    )

    # work out beam waist here from HG decomposition

    x = np.linspace(0, nx * pixel_size, 500)
    X, Y = np.meshgrid(x, x)
    transverse_field = transverse_profile.evaluate(X,Y)

    # evaluate transverse profile to get field to pass into this

    waist = estimate_best_HG_waist(x, x, transverse_field)

    if dropout_figs:

        # Re-evaluate transverse profile with figure extend determined from waist from hermite-gauss decomposition
        pltextent = np.array([-5*waist, 5*waist, -5*waist, 5*waist])   
        x = np.linspace(-5*waist, 5*waist, 500)
        X, Y = np.meshgrid(x, x)
        transverse_profile_evaluated = np.abs(transverse_profile.evaluate(X,Y))**2

        plt.imshow(transverse_profile_evaluated, cmap="magma", extent = pltextent)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("LASY transverse_profile_from_data")
        plt.show()

    return transverse_profile, waist

def densityplot3d(ax: plt.axis, 
                  x: np.ndarray, 
                  y: np.ndarray, 
                  z: np.ndarray, 
                  values: np.ndarray, 
                  decay: float = 2.0,
                  opacity: float = 0.7, 
                  cmap: mpl.colors.LinearSegmentedColormap = plt.cm.jet,
                  **kwargs) -> None:
    """
    Create a density plot for X, Y, Z coordinates and corresponding intensity values.

    Parameters:
    -----------
    ax : plt.axis
        The axis object to plot the density plot on.
    x : np.ndarray
        An array of X-coordinates.
    y : np.ndarray
        An array of Y-coordinates.
    z : np.ndarray
        An array of Z-coordinates.
    values : np.ndarray
        An array of intensity values.
    decay : float, optional
        The decay factor for the alpha values. Default is 2.0.
    opacity : float, optional
        The opacity value for the alpha values. Default is 1.0.
    cmap : mpl.colors.LinearSegmentedColormap, optional
        The colormap used for mapping intensity values to RGB colors. Default is plt.cm.jet.
    **kwargs
        Additional keyword arguments to pass to the scatter function.

    Returns:
    --------
    None
    """

    # Calculate RGB colors from intensities
    # Normalize the intensities between 0 and 1 and convert them to RGB colors using the chosen colormap
    normed_values = values / np.max(values)
    colors = cmap(normed_values)
    
    # Create alpha values for each data point based on its intensity and the specified decay factor
    alphas = (values / np.max(values)) ** decay
    alphas *= opacity
    
    colors[:, :, :, 3] = alphas  # add alpha values to RGB values
    
    # Flatten color array but keep last dimension
    colors_flattened = colors.reshape(-1, colors.shape[-1])  

    # Plot a 3D scatter with adjusted alphas
    ax.scatter(x, y, z, c=colors_flattened, **kwargs)
    
    return None


def ProfileFromInsight(insight_file, polarisation, wavelength, dropout_figs):

    laser_profile = FromInsightFile(insight_file, polarisation, omega0=wavelength) # omega0 = barycenter by default

    if dropout_figs:

        # FIND BETTER WAY TO DEFINE THESE RANGES BASED ON LASER PULSE
        xy_range = (-0.5e-4, 0.5e-4)
        z_range = (-40e-15, 40e-15)

        # create grid to evaluate laser_profile on
        x_array = np.linspace(-50e-6, 50e-6, 100)
        y_array = np.linspace(-50e-6, 50e-6, 100)
        t_array = np.linspace(-40e-15, 40e-15, 100)

        X,Y,T = np.meshgrid(x_array, y_array, t_array)

        # CHANGE DENSITYPLOT3DFUNCTION TO SET THESE AUTOMATICALLY
        # Set up parameters for plot display
        opacity = 0.6
        decay = 3
        labels = ['x', 'y', 't']
        cmap = plt.cm.hot
        markersize = 3

        # Create figure and axis
        fig = plt.figure(num=1, clear=True, figsize=(12, 6))
        ax = fig.add_subplot(111, projection='3d')

        laser_profile_evaluated = laser_profile.evaluate(X,Y,T)

        laser_profile_amplitude = np.abs(laser_profile_evaluated)**2
        # SHOULD i SQUARE AMPLITUDE OR NOT?

        # Plot data
        densityplot3d(ax, X, T, Y, laser_profile_amplitude, cmap=cmap, opacity=opacity, decay=decay, s=markersize)

        # Set up colorbar
        cmap = plt.cm.ScalarMappable(norm=mpl.colors.Normalize(np.min(laser_profile_amplitude), np.max(laser_profile_amplitude)), cmap=cmap)

        # Add colorbar with labels and padding
        fig.colorbar(cmap, ax=ax, label='Intensity', pad=0.1)

        # Set labels and add padding
        ax.set_xlabel(labels[0], labelpad=5)
        ax.set_zlabel(labels[1], labelpad=5)
        ax.set_ylabel(labels[2], labelpad=5)

        # Set limits
        adjuster = 1.05
        ax.set_xlim([xy_range[0] * adjuster, xy_range[1] * adjuster])
        ax.set_zlim([xy_range[0] * adjuster, xy_range[1] * adjuster])
        ax.set_ylim([z_range[0] * adjuster, z_range[1] * adjuster])

        # Adjust tick frequency
        stepsize = 0.5
        ax.xaxis.set_ticks(np.linspace(xy_range[0], xy_range[1], 5))
        ax.zaxis.set_ticks(np.linspace(xy_range[0], xy_range[1], 5))
        ax.yaxis.set_ticks(np.linspace(z_range[0], z_range[1], 5))

        ax.set_title("Laser pulse from INSIGHT measurement")

        # Optional change of view port
        # ax.view_init(30, 120)

        # Save plot
        #fig.savefig('Density_plot.png', bbox_inches='tight', dpi=300)

        plt.show()


    return laser_profile