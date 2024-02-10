import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import logging

# Set up logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger()

# Get fits data
# main dimension (2602, 3906)
def get_fits_blue_data(folder_name):
    directory = f"../data/{folder_name}/"
    images_blue_data = []
    wrong_files_counter = 0
    if os.path.exists(directory):
        logger.info('Getting blue images data...')
        for file_name in os.listdir(directory):
            temp_list = fits.open(directory + file_name)
            temp_image = np.dot(temp_list[0].data.T, [1, 1, 1])
            if temp_image.shape == (3906, 2602):
                images_blue_data.append(temp_image)
            else:
                logger.info(f'{file_name} is not (2602, 3906) shape...')
                wrong_files_counter += 1
            temp_list.close()
        logger.info(f'{wrong_files_counter} images were not in the correct format.')
        logger.info('# get_fits_blue_data # worked successfully...')
        return images_blue_data
    else:
        print("File does not exist.")

# Master Dark process:
def get_master_dark(folder_name):
    logger.info('Master dark processing started...')
    images_blue_data = get_fits_blue_data(folder_name)
    staked_blue_images = np.stack(images_blue_data, axis=0)
    median_array = np.median(staked_blue_images, axis=0)
    std_array = np.std(staked_blue_images, axis=0)
    master_dark_image_shape = staked_blue_images.shape[1], staked_blue_images.shape[2]
    master_dark_image = np.zeros(master_dark_image_shape)
    for i in range(staked_blue_images.shape[1]):
        for j in range(staked_blue_images.shape[2]):
            temp_median = median_array[i][j]
            temp_std = std_array[i][j]
            pixels_value = []
            for image in images_blue_data:
                if temp_median - 3 * temp_std < image[i][j] < temp_median + 3 * temp_std:
                    pixels_value.append(image[i][j])
            if pixels_value:
                master_dark_image[i][j] = np.mean(np.array(pixels_value))
    logger.info('Master dark processing worked successfully...')
    return master_dark_image

# Master flat process:
def get_master_flat(folder_name):
    logger.info('Master flat processing started...')
    images_blue_data = get_fits_blue_data(folder_name)
    staked_blue_images = np.stack(images_blue_data, axis=0)
    median_array = np.median(staked_blue_images, axis=0)
    std_array = np.std(staked_blue_images, axis=0)
    master_flat_image_shape = staked_blue_images.shape[1], staked_blue_images.shape[2]
    master_flat_image = np.zeros(master_flat_image_shape)
    for i in range(staked_blue_images.shape[1]):
        for j in range(staked_blue_images.shape[2]):
            temp_median = median_array[i][j]
            temp_std = std_array[i][j]
            pixels_value = []
            for image in images_blue_data:
                if temp_median - 3 * temp_std < image[i][j] < temp_median + 3 * temp_std:
                    pixels_value.append(image[i][j])
            if pixels_value:
                master_flat_image[i][j] = np.mean(np.array(pixels_value))
    logger.info('Getting master flat dark...')
    master_flat_dark = np.load("../data/master_darks/master_dark_flat_evening_fits.npy")
    master_flat_image = master_flat_image - master_flat_dark
    logger.info('Master flat processing worked successfully...')
    return master_flat_image

# Gain table process:
def get_gain_table(folder_name):
    logger.info('Gain table processing started...')
    master_flat_image = get_master_flat(folder_name)
    gain_table = master_flat_image / np.median(master_flat_image)
    logger.info('Gain table processing worked successfully...')
    return gain_table    
    
# Process multiple folders for master dark images
folder_names =[
    "dark_flat_evening_fits",
    "dark_30s_iso200_fits",
    "dark_30s_iso1600_fits"
]
for folder_name in folder_names:
    master_dark_image = get_master_dark(folder_name)
    np.save(f"../data/master_darks/master_{folder_name}.npy", master_dark_image)
    logger.info(f"master_{folder_name}.npy saved successfully...\n\n\n")
    
# Process a single folder for gain table
folder_name = "flat_evening_fits"
gain_table = get_gain_table(folder_name)
np.save("../data/gain_table/gain_table.npy", gain_table)
logger.info("gain_table.npy saved successfully...")
