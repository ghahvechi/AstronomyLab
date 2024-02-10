import logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.optimize import curve_fit

# Set up logging
logging.basicConfig(level=logging.NOTSET)
logger = logging.getLogger()

# Opening FITS file
logger.info('Opening FITS file...')
image = fits.open("D:\Files\PSF\stars.fits")[0].data
logger.info('Opened successfully')

# Identify stars
number_of_stars = 0
stars_0 = []
stars_1 = []
median = np.median(image)
std = np.std(image)
for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        if image[i][j] > 60000 and 400 < i < 1556 and 400 < j < 2321:  
            check = 1
            if number_of_stars != 0:
                for c in range(number_of_stars):
                    if abs(i - stars_0[c]) < 40 and abs(j - stars_1[c]) < 40:
                        check = 0  
            if check == 1:  
                stars_0.append(i)
                stars_1.append(j)
                number_of_stars += 1

# Define functions for center of stars and standard deviation calculation
def center_of_stars(stars_0, stars_1, number_of_stars, image):
    # Implementation omitted for brevity
    pass

def process_every_crop(number_of_stars, i_centers, j_centers, maxes, image):
    # Implementation omitted for brevity
    pass

def gaussian_function(x, a, x0, sigma, b, c):
    # Implementation omitted for brevity
    pass

def fit(values_of_crope, number_of_stars):
    # Implementation omitted for brevity
    pass

# Calculate standard deviation for each star
stds = []
for i in range(number_of_stars):
    stds.append(fit(process_every_crop(number_of_stars, center_of_stars(stars_0, stars_1, number_of_stars, image)[0], center_of_stars(stars_0, stars_1, number_of_stars, image)[1], center_of_stars(stars_0, stars_1, number_of_stars, image)[2], image)[i], number_of_stars))

# Plot PSF of all stars in one
# Plot histogram of standard deviations
# Plot PSF of each star individually

# Reporting standard deviation
my_std = np.std(stds)
median = np.median(stds)

new_stds = []

for std in stds:
    if std <= median + 2 * my_std:
        new_stds.append(std)
print(np.mean(new_stds))
print(np.std(new_stds))

# Print table
for i in range(number_of_stars):
    print(center_of_stars(stars_0, stars_1, number_of_stars, image)[0][i], center_of_stars(stars_0, stars_1, number_of_stars, image)[1][i], center_of_stars(stars_0, stars_1, number_of_stars, image)[2][i], stds[i])
