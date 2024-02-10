import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

# Function to convert 3D data to 2D
def convert_3d_to_2d(three_D_img_datas):
    images = []
    for three_D_img_data in three_D_img_datas:
        # Convert 3D data to 2D using dot product with color weights
        temp_2D_data = np.dot(three_D_img_data.T, [0.2989, 0.5870, 0.1140])
        images.append(temp_2D_data)
    return images

# Function to extract data from FITS files
def get_fits_data(folder_name):
    directory = f"D:\\Files\\{folder_name}\\"
    three_D_img_datas = []
    for filename in os.listdir(directory):
        temp_list = fits.open(directory + filename)
        three_D_img_datas.append(temp_list[0].data)
        temp_list.close()
    return three_D_img_datas

# Function to get the dark image
def get_dark_image(folder_name):
    # Convert 3D data to 2D images
    images = convert_3d_to_2d(get_fits_data(folder_name))

    # Stack images along the depth axis
    staked_images = np.stack(images, axis=0)

    # Compute median and standard deviation arrays
    median_array = np.median(staked_images, axis=0)
    std_array = np.std(staked_images, axis=0)

    # Initialize dark image array
    dark_image_shape = staked_images.shape[1], staked_images.shape[2]
    dark_image = np.zeros(dark_image_shape)

    # Iterate over pixels
    for i in range(staked_images.shape[1]):
        for j in range(staked_images.shape[2]):
            temp_median = median_array[i][j]
            temp_std = std_array[i][j]
            pixels_value = []
            # Select pixels within 3 standard deviations of the median
            for image in images:
                if temp_median - 3 * temp_std < image[i][j] < temp_median + 3 * temp_std:
                    pixels_value.append(image[i][j])
            # Compute median of selected pixels and assign to dark image
            if pixels_value:
                dark_image[i][j] = np.median(np.array(pixels_value))
    return dark_image
