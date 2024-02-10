import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

# Configure logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger()

# Function to extract data from FITS files
def get_fits_data(folder_name):
    logger.info('Converting to FITS file...')
    directory = f"D:\\Files\\dark_set2_fits\\{folder_name}\\" # Target location
    three_D_img_datas = []
    for filename in os.listdir(directory):
        temp_list = fits.open(directory + filename)
        three_D_img_datas.append(temp_list[0].data)
        temp_list.close()
    logger.info('Conversion successful')
    return three_D_img_datas

# Function to convert 3D data to 2D
def convert_3d_to_2d(three_D_img_datas):
    logger.info('Image processing started...')
    images = []
    for three_D_img_data in three_D_img_datas:
        temp_2D_data = np.dot(three_D_img_data.T, [0.2989, 0.5870, 0.1140])
        images.append(temp_2D_data)
    return images

# Function to get dark image, median, and standard deviation
def get_dark_image(folder_name):
    images = convert_3d_to_2d(get_fits_data(folder_name))

    logger.info('Stacking images...')
    staked_images = np.stack(images, axis=0)
    median_array = np.median(staked_images, axis=0)
    std_array = np.std(staked_images, axis=0)

    dark_image_shape = staked_images.shape[1], staked_images.shape[2]
    dark_image = np.zeros(dark_image_shape)
    values_in_std_range = []
    for i in range(staked_images.shape[1]):
        for j in range(staked_images.shape[2]):
            temp_median = median_array[i][j]
            temp_std = std_array[i][j]
            pixels_value = []
            for image in images:
                if temp_median - 3 * temp_std < image[i][j] < temp_median + 3 * temp_std:
                    pixels_value.append(image[i][j])
                    values_in_std_range.append(image[i][j])
            if pixels_value:
                dark_image[i][j] = np.median(np.array(pixels_value))
    values_in_std_range_array = np.array(values_in_std_range)
    median = np.median(values_in_std_range_array)
    std = np.std(values_in_std_range_array)
    return dark_image, median, std

# Function to draw histogram
def draw_histogram(dark_image, exposure_time):
    logger.info('Drawing histogram...')
    title = exposure_time + ' s'
    file_path = 'D:\\Files\\Dark\\Histogram\\' + exposure_time.replace('/', '_') + '.png' # Target location
    bins_list = list(map(lambda x: x * 0.5, range(7 * 2)))
    counts, bins = np.histogram(dark_image, bins=bins_list)
    plt.clf()
    plt.hist(bins[:-1], bins, weights=counts)
    plt.title(title)
    plt.xlabel("Value of pixel")
    plt.ylabel("Number of pixels")
    plt.savefig(file_path)

# Dictionary of exposure times and corresponding folder names
folders_name_dic = {"1/4000": "1_over_4000", "1/1000": "1_over_1000", "1/200": "1_over_200",
                    "1/100": "1_over_100", "1/30": "1_over_30",
                    "1/10": "1_over_10", "30": "30", "10": "10", "1": "1"}

# Lists to store medians and standard deviations
medians_list = []
stds_list = []

# Iterate over exposure times
for exposure_time in folders_name_dic.keys():
    logger.info(f'Exposure time = {exposure_time}')
    dark_image, median, std = get_dark_image(folders_name_dic[exposure_time])
    draw_histogram(dark_image, exposure_time)
    medians_list.append(median)
    stds_list.append(std)
    logger.info(f"{exposure_time} : median = {median} || std = {std}")
    # Write results to a text file
    with open("D:\\Files\\Dark\\output.txt", 'a') as file:
        file.write(f"{exposure_time} : median = {median} || std = {std}" + "\n")

# Draw medians plot
xs = np.log10([1 / 4000, 1 / 1000, 1 / 200, 1 / 100, 1 / 30, 1 / 10, 30, 10, 1])
medians_list = np.array([3.9072, 3.9287, 3.8885, 3.8857, 3.9287, 3.847, 3.3741, 3.6174, 3.8947])
plt.clf()
plt.scatter(xs, medians_list)
plt.title('Signal Medians')
plt.ylabel('Signal Median')
plt.xlabel('Log(Time Delta)')
plt.savefig('D:\\Files\\Dark\\medians.png') # Save in the target Location

# Draw standard deviations plot
stds_list = np.array([4.452546064522884, 4.431540601977413, 4.407476891371636, 4.412552303580464, 4.44396557578929,
                      4.439988864992894, 19.48276329993115, 13.865098753867702, 7.499661071179732])
plt.clf()
plt.scatter(xs, stds_list)
plt.title('Signal STDs')
plt.ylabel('Signal STD')
plt.xlabel('Log(Time Delta)')
plt.savefig('D:\\Files\\Dark\\stds.png') # Save in the target location
