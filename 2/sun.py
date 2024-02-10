import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import dark_proccesing as dp  # Importing a module named "dark_proccesing" as "dp"
from scipy.stats import linregress

def sun_data_dark_current_correction(sun_image_uncorrected, dark_image):
    """
    Subtract dark current from sun image data.

    :param sun_image_uncorrected: ndarray
    :param dark_image: ndarray
    :return: Corrected sun image data
    """
    if dark_image.shape == sun_image_uncorrected.shape:
        # Initialize array for corrected sun image
        sun_image_corrected_shape = sun_image_uncorrected.shape[0], sun_image_uncorrected.shape[1]
        sun_image_corrected = np.zeros(sun_image_corrected_shape)

        # Subtract dark data from sun data pixel-wise
        for i in range(dark_image.shape[0]):
            for j in range(dark_image.shape[1]):
                sun_image_corrected[i][j] = sun_image_uncorrected[i][j] - dark_image[i][j]
        return sun_image_corrected
    else:
        print("!!! Dark image shape and sun image shape are not the same !!!")

def get_center_of_mass(image):
    """
    Compute the position of the center of mass of the given image.

    :param image: ndarray representing the sun image
    :return: Array containing the position of the center of mass
    """
    total_mass = 0
    weight_vector = np.array([0., 0.])

    # Iterate over pixels to calculate total mass and weighted sum
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            weight_vector += image[i][j] * np.array([i, j])
            total_mass += image[i][j]

    # Calculate the center of mass
    center_of_mass = np.round(weight_vector / total_mass)
    return center_of_mass

def get_fwhm(data):
    """
    Compute the Full Width at Half Maximum (FWHM) of the given data.

    :param data: ndarray representing the horizontal data of the sun image
    :return: FWHM value
    """
    max_data = np.max(data)
    half_max_data = max_data / 2

    # Find indices where data is greater than or equal to half maximum
    indices = np.where(data >= half_max_data)[0]

    # Compute FWHM
    fwhm = (indices.max() - indices.min())
    return fwhm

def sun_radial_profile(data, i_cm):
    """
    Generate and save the radial profile plot of the sun image.

    :param data: ndarray representing the horizontal data of the sun image
    :param i_cm: Y-coordinate of the center of mass
    """
    title = "Sun Profile"
    x_label = "Position"
    y_label = "Intensity"

    # Compute relative positions
    x_data = np.where(data)[0] - i_cm
    normal_data = data / np.max(data)

    # Plot and save radial profile
    plt.clf()
    plt.plot(x_data, normal_data, '-')
    plt.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig('D:\\Files\\Sun\\sun_radial_profile.png')

def eddington_plot(data, i_cm, radius):
    """
    Generate and save the Eddington plot of the sun image.

    :param data: ndarray representing the horizontal data of the sun image
    :param i_cm: Y-coordinate of the center of mass
    :param radius: Radius of the sun
    :return: R-squared value of the linear regression
    """
    x_label = "μ = cos(θ)"
    y_label = "Intensity"

    # Select sun data within the specified radius
    sun_data = data[i_cm - radius:i_cm + 1] / np.max(data)
    mu = np.cos(np.linspace(np.deg2rad(90), 0, num=sun_data.shape[0]))

    # Compute Eddington relation
    I_mu = (data[i_cm] / np.max(data)) * (2 + 3 * mu) / 5

    # Plot sun radial intensity and Eddington relation
    plt.clf()
    plt.plot(mu, sun_data, '-', label="Sun Radial Intensity")
    plt.plot(mu, I_mu, 'r-', label="Eddington Relation")
    plt.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.savefig('D:\\Files\\Sun\\eddington_plot.png')

    # Perform linear regression and return R-squared value
    slope, intercept, r_value, p_value, std_err = linregress(I_mu, sun_data)
    return r_value

def eddington_plot_2(data, i_cm, radius):
    """
    Generate and save the second type of Eddington plot of the sun image.

    :param data: ndarray representing the horizontal data of the sun image
    :param i_cm: Y-coordinate of the center of mass
    :param radius: Radius of the sun
    """
    x_label = "θ"
    y_label = "Intensity"

    # Select sun data within the specified radius
    sun_data = data[i_cm - radius:i_cm + 1] / np.max(data)
    theta = np.linspace(90, 0, num=sun_data.shape[0])
    mu = np.cos(np.linspace(np.deg2rad(90), 0, num=sun_data.shape[0]))
    I_mu = (data[i_cm] / np.max(data)) * (2 + 3 * mu) / 5

    # Plot sun radial intensity and Eddington relation
    plt.clf()
    plt.plot(theta, sun_data, '-', label="Sun Radial Intensity")
    plt.plot(theta, I_mu, 'r-', label="Eddington Relation")
    plt.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.savefig('D:\\Files\\Sun\\eddington_plot_2.png')

def radial_profile_outside_sun(data, i_cm, radius):
    """
    Generate and save the radial profile outside the sun.

    :param data: ndarray representing the horizontal data of the sun image
    :param i_cm: Y-coordinate of the center of mass
    :param radius: Radius of the sun
    """
    title = "Outside Radial Profile"
    x_label = "Position"
    y_label = "Intensity"

    # Select data outside the sun
    outside_data = data[:i_cm - radius] / np.max(data)
    x_data = np.where(outside_data)[0] - i_cm

    # Plot and save radial profile outside the sun
    plt.clf()
    plt.plot(x_data, outside_data, '-')
    plt.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig('D:\\Files\\Sun\\radial_profile_outside_sun.png')

# Load corrected sun image data
sun_image_corrected = np.load('sun_corrected.npy')

# Calculate center of mass and sun radius
location = get_center_of_mass(sun_image_corrected)
horizontal_data = sun_image_corrected[int(location[0])]
sun_radius = int(np.round(get_fwhm(horizontal_data)) / 2) + 1
y_cm = int(location[1])

# Generate and save plots
sun_radial_profile(horizontal_data, y_cm)
r_value = eddington_plot(horizontal_data, y_cm, sun_radius)
eddington_plot_2(horizontal_data, y_cm, sun_radius)
radial_profile_outside_sun(horizontal_data, y_cm, sun_radius)

# Write results to output file
with open("D:\\Files\\Sun\\output.txt", 'w') as file:
    file.write(f"Location: ({location[0]}, {y_cm})\nSun Radius: {sun_radius}\nR-squared: {r_value}")
