# AstronomyLab

## Converting Raw Data to FITS
To convert raw images from a camera to FITS format, we utilize the dcraw and convert commands. Below is a simple script that automates this conversion process:
### Converter Script
```bash
raw_dir=raw
fits_dir=fits

# Create FITS directory if it doesn't exist
mkdir -p $fits_dir

for raw_filename in "$raw_dir"/*; do
  # Determine output FITS filename
  fits_filename="${raw_filename%.NEF}.FITS"
  fits_filename="${fits_filename%.CR2}.FITS"
  fits_filename="${fits_filename/raw/fits}"

  # Convert raw to FITS
  dcraw -c -w -v -4 -S 65535 $raw_filename | pamtofits > $fits_filename
done

```
### Usage
1. Place your raw image files in the raw directory.
2. Run the script by executing it in your terminal.
3. Converted FITS files will be saved in the fits directory.

### Notes
. This script assumes that your raw files have either a .NEF or .CR2 extension. Adjust the script if your raw files have a different extension.
. Make sure to have the dcraw and pamtofits utilities installed on your system for this script to work properly.
. You can customize the script further based on your specific requirements or preferences.
