#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit
from PIL import Image

# Setting up logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger()

# Function to get the center of mass of a star
def get_center_of_mass(image, i_of_star, j_of_star, radius):
    """
    Calculate the position of the center of mass.

    :param image: ndarray representing the image
    :param i_of_star: row index of the star
    :param j_of_star: column index of the star
    :param radius: radius around the star to consider
    :return: 2-element array representing the position of the center of mass
    """
    total_mass = 0
    weight_vector = np.array([0., 0.])
    for i in range(i_of_star - radius, i_of_star + radius + 1):
        for j in range(j_of_star - radius, j_of_star + radius + 1):
            weight_vector += image[i][j] * np.array([i, j])
            total_mass += image[i][j]
    center_of_mass = (np.round(weight_vector / total_mass)).astype(int)
    return center_of_mass

# Function to calculate the signal of a star
def calculate_signal(image, i_of_star, j_of_star, radius, estimated_radius):
    signal = signal_minus_sky = pixel_number = sky = 0
    for i in range(i_of_star - radius, i_of_star + radius + 1):
        for j in range(j_of_star - radius, j_of_star + radius + 1):
            x_val = i - i_of_star
            y_val = j - j_of_star
            if (x_val ** 2 + y_val ** 2 <= radius ** 2):
                signal += image[i][j]
                pixel_number += 1
    sky = sky_mean(image, i_of_star, j_of_star, estimated_radius)
    signal_minus_sky = signal - (pixel_number * sky)
    return signal, signal_minus_sky, pixel_number

# Function to calculate the mean sky background
def sky_mean(image, i_of_star, j_of_star, estimated_radius):
    noise = 0
    pixel_counter = 0
    big_radius = np.round(estimated_radius * 2.5).astype(int)
    for i in range(i_of_star - big_radius, i_of_star + big_radius + 1):
        for j in range(j_of_star - big_radius, j_of_star + big_radius + 1):
            x_val = i - i_of_star
            y_val = j - j_of_star
            if ((x_val ** 2 + y_val ** 2 >= (2 * estimated_radius) ** 2) and (x_val ** 2 + y_val ** 2 <= (2.5 * estimated_radius) ** 2)):
                noise += image[i][j]
                pixel_counter += 1
    return noise / pixel_counter

# Function to calculate the radius of a star
def calculate_radius_of_star(image, i_of_star, j_of_star, estimated_radius):
    snr = []
    signal = noise = SNR = max_SNR = radius = pixel_number = sky = signal_minus_sky = best_signal_minus_sky = 0
    for i in range(1, 30):
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

# Function to plot SNR
def plot_SNR(SNR, name_of_star):
    x = list(range(29))
    title = name_of_star
    plt.clf()
    plt.plot(x, SNR, 'b-')
    SNR_max = np.max(SNR)
    max_index = np.where(SNR == SNR_max)[0][0]
    plt.plot([x[max_index], x[max_index]], [0, SNR_max], color='r', linestyle=':')
    plt.grid(color='r', alpha=0.5, linestyle='dashed', linewidth=0.5)
    plt.title(title)
    plt.ylabel('SNR')
    plt.xlabel('radius')
    plt.savefig(f'D:\\Asnave\\evening\\M29_SNRs\\{title}.png')

# Function to find the centers of stars
def center_of_stars(i_stars, j_stars, number_of_stars, image):
    i_centers = []
    j_centers = []
    maxes = []
    for counter in range(number_of_stars):
        temp_max = image[i_stars[counter]][j_stars[counter]]
        i = (i_stars[counter] - 25)
        j = (j_stars[counter] - 25)
        temp_i = i_stars[counter]
        temp_j = j_stars[counter]
        while i < (i_stars[counter] + 25):
            while j < (j_stars[counter] + 25):
                if image[i][j] > temp_max:
                    temp_max = image[i][j]
                    temp_i = i
                    temp_j = j
                j += 1
            i += 1
            j = (j_stars[counter] - 40)
        maxes.append(temp_max)
        i_centers.append(temp_i)
        j_centers.append(temp_j)
    return i_centers, j_centers, maxes

# Load image data
image = np.load("D:\\Asnave\\evening\\final_stars_image\\final_stars_image.npy")

# Identifying the stars
number_of_stars = 0
i_stars = []
j_stars = []
median = np.median(image)
std = np.std(image)
for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        if 3500 < image[i][j] and 60 < i < 1891 and 60 < j < 1387:  # Adjusted condition to filter stars
            check = 1
            if number_of_stars != 0:
                for c in range(number_of_stars):
                    if abs(i - i_stars[c]) < 40 and abs(j - j_stars[c]) < 40:
                        check = 0  # Avoid counting the same star twice
            if check == 1:
                i_stars.append(i)
                j_stars.append(j)
                number_of_stars += 1

# Refine the centers of the stars
i_centers, j_centers, maxes = center_of_stars(i_stars, j_stars, number_of_stars, image)
for k in range(number_of_stars):
    COM = get_center_of_mass(image, i_centers[k], j_centers[k], 17)
    i_centers[k] = COM[0]
    j_centers[k] = COM[1]

# Calculations for each star
i_base_COM = 248  # Base star position
j_base_COM = 988
base_max_SNR, base_SNR, base_radius, base_signal_minus_sky = calculate_radius_of_star(image, i_base_COM, j_base_COM, 15)
base_error = 1.0875 / base_max_SNR
base_magnitude = 9.94
magnitudes = []
radii = []
SNRs = []
errors = []
last_i = []
last_j = []
counter = 1
for i in range(number_of_stars):
    temp_SNR, temp_SNRs, temp_radius, temp_signal_minus_sky = calculate_radius_of_star(image, i_centers[i], j_centers[i], 15)
    if 5 < temp_radius < 18:  # Adjusted condition to filter stars by radius
        plot_SNR(temp_SNRs, counter)
        temp_error = 1.0875 / temp_SNR
        error = np.sqrt((base_error ** 2) + (temp_error ** 2))
        SNRs.append(temp_SNR)
        radii.append(temp_radius)
        last_i.append(i_centers[i])
        errors.append(error)
        temp_magnitude = -2.5 * np.log10(temp_signal_minus_sky / base_signal_minus_sky) + base_magnitude
        magnitudes.append(temp_magnitude)
        counter += 1
