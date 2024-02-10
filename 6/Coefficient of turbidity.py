#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries
import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
from PIL import Image

# Set up logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger()

# Reading FITS files
logger.info('Getting FITS files of base star, dark, and gain table...')
three_D_image_of_base = fits.open("D:\\Asnave\\evening\\stars_evening_fits\\base_onogh\\IMG_8906_14020427_10_13_45.FITS")[0].data         
master_dark_iso200 = np.load("D:\\Asnave\\evening\\master_darks\\master_dark_30s_iso200_fits.npy")
gain_table = np.load("D:\\Asnave\\evening\\gain_table\\gain_table.npy")
logger.info('Files opened successfully.')

# Converting 3D to 2D for base, dark, and flat processing
logger.info('Converting 3D to 2D for base, dark, and flat processing...')
base_image = np.dot(three_D_image_of_base.T, [1, 1, 1])
base_image = (base_image - master_dark_iso200) / gain_table
logger.info('Base image is 2D and dark and flat processing analyzed successfully.')

# Get images and process them
images = []
wrong_file_counter = 0
number_of_images = 0
folder = 'D:\\Asnave\\evening\\stars_evening_fits\\Onogh\\'
if os.path.exists(folder):
    logger.info('Getting FITS files...')
    for file_name in os.listdir(folder):
        temp_list = fits.open(folder + file_name)[0].data 
        temp_image = np.dot(temp_list.T, [1, 1, 1])
        temp_image = (temp_image - master_dark_iso200) / gain_table
        if temp_image.shape == (3906, 2602):
            images.append(temp_image)
            logger.info(f'{file_name} added.')
            number_of_images += 1
        else:
            logger.info(f'{file_name} is not (3906, 2602) shape.')
            wrong_file_counter += 1
    logger.info(f'{wrong_file_counter} images were not in correct format.')
    logger.info('Getting FITS data worked successfully.')
else:
    logger.info('Folder does not exist.')

# Coordinates of the base star and other stars
i_base_COM = 1662
j_base_COM = 759
i_COM = [1869, 1864, 1839, 1883, 2127, 2123, 2071, 2043, 2027, 2002, 1897, 1881, 1832, 1832, 1732, 1728]                
j_COM = [681, 706, 741, 717, 616, 643, 680, 604, 648, 637, 608, 582, 610, 630, 1534, 1502]

# Function to calculate signal within a radius
def calculate_signal(image, i_of_star, j_of_star, radius, estimated_radius):
    signal = signal_minus_sky = pixel_number = sky = 0
    for i in range(i_of_star - radius, i_of_star + radius + 1):
        for j in range(j_of_star - radius, j_of_star + radius + 1):
            x_val = i - i_of_star
            y_val = j - j_of_star
            if (x_val**2 + y_val**2 <= radius**2):
                signal += image[i][j]
                pixel_number += 1
    sky = sky_mean(image, i_of_star, j_of_star, estimated_radius)
    signal_minus_sky = signal - (pixel_number * sky)
    return signal, signal_minus_sky, pixel_number

# Function to calculate the mean sky value
def sky_mean(image, i_of_star, j_of_star, estimated_radius):
    noise = 0
    pixel_counter = 0
    big_radius = np.round(estimated_radius * 2.5).astype(int)
    for i in range(i_of_star - big_radius, i_of_star + big_radius + 1):
        for j in range(j_of_star - big_radius, j_of_star + big_radius + 1):
            x_val = i - i_of_star
            y_val = j - j_of_star
            if ((x_val**2 + y_val**2 >= (2 * estimated_radius)**2) and (x_val**2 + y_val**2 <= (2.5 * estimated_radius)**2)):
                noise += image[i][j]
                pixel_counter += 1
    return noise / pixel_counter

# Function to get the center of mass of a star
def get_center_of_mass(image, i_of_star, j_of_star):
    total_mass = 0
    weight_vector = np.array([0., 0.])
    for i in range(i_of_star - 20, i_of_star + 21):
        for j in range(j_of_star - 20, j_of_star + 21):
            weight_vector += image[i][j] * np.array([i, j])
            total_mass += image[i][j]
    center_of_mass = (np.round(weight_vector / total_mass)).astype(int)
    return center_of_mass

# Function to plot SNR
def plot_SNR(SNR):
    x = list(range(100))
    title = "18 SNR for Unukalhai at 12dsfeawfdasfas_03_15"
    plt.clf()
    plt.plot(x, SNR, 'b-')
    SNR_max = np.max(SNR)
    max_index = np.where(SNR == SNR_max)[0][0]
    plt.plot([x[max_index], x[max_index]], [0, SNR_max], color='r', linestyle=':')
    plt.grid(color='r', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.title(title)
    plt.ylabel('SNR')
    plt.xlabel('radius')
    plt.savefig(f'D:\\Asnave\\evening\\Onogh_SNR\\{title}.png')

# Function to calculate the radius of a star
def calculate_radius_of_star(image, i_of_star, j_of_star, estimated_radius):
    snr = []
    signal = noise = SNR = max_SNR = radius = pixel_number = sky = signal_minus_sky = best_signal_minus_sky = 0
    for i in range(1, 101):
        signal, signal_minus_sky, pixel_number = calculate_signal(image, i_of_star, j_of_star, i, estimated_radius)            
        temp_noise = pixel_number + signal
        noise = np.sqrt(temp_noise)
        SNR = signal_minus_sky / noise
        snr.append(SNR)
        if (SNR > max_SNR):
            max_SNR = SNR
            radius = i
            best_signal_minus_sky = signal_minus_sky
    return max_SNR, snr, radius, best_signal_minus_sky

# Calculate base star parameters
base_max_SNR, base_SNR, base_radius, base_signal_minus_sky = calculate_radius_of_star(base_image, i_base_COM, j_base_COM, 15)                  
base_error = 1.0875 / base_max_SNR
base_magnitude = 0
magnitudes = []
radii = []
SNRs = []
error = []
signal_minus_skys = []
errors = []

# Calculate parameters for other stars
for i in range(number_of_images):
    temp_SNR, temp_SNRs, temp_radius, temp_signal_minus_sky = calculate_radius_of_star(images[i], i_COM[i], j_COM[i], 15)                  
    temp_error = 1.0875 / temp_SNR
    error = np.sqrt((base_error**2) + (temp_error**2))
    SNRs.append(temp_SNR)
    radii.append(temp_radius)
    signal_minus_skys.append(temp_signal_minus_sky)
    errors.append(error)
    temp_magnitude = -2.5 * np.log10(temp_signal_minus_sky / base_signal_minus_sky) + base_magnitude
    magnitudes.append(temp_magnitude)

# Calculate error of slope
sec_z = np.array([1.362899411, 1.399633644, 1.447424314, 1.533510884, 1.555723827, 2.069183638, 2.082363744])
magnitude = np.array([-0.015056194, -0.003240927, 0.003984195, 0.064878515, 0.073275622, 0.317885486, 0.369352119])
A = np.vstack([sec_z, np.ones(len(sec_z))]).T
m, b = np.linalg.lstsq(A, magnitude, rcond=None)[0]

# Plot fitted line
plt.plot(sec_z, magnitude, "ob")
plt.plot(sec_z, m * sec_z + b, 'r', label="Fitted line")
plt.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)
plt.xlabel("Relative magnitude")
plt.ylabel("Sec (z)")
plt.legend()
plt.savefig(f'D:\\Asnave\\fitted_line.png')
print(f"m = {m} and b = {b}")

# Calculate error of slope
x = sec_z
y = magnitude
d = y - m * x - b
x_bar = np.sum(x) / len(x)
D = np.sum(np.square(x - x_bar))
delta_a = np.sqrt((1 / D) * (np.sum(np.square(d)) / (len(x) - 2)))
delta_b = np.sqrt(((1 / len(x)) + (np.square(x_bar) / D)) * (np.sum(np.square(d)) / (len(x) - 2)))

print(f"dleta m = {delta_a} and delta_b = {delta_b}")
