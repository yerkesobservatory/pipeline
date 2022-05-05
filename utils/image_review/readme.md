# Installing
Just run `pip install -r requirements.txt` to install the necessary Python dependencies.

# Running
To run the image review script
```
python3 image_review.py --outfolder <folder to save your selected images to> <input files>
```

For example, if I had a list of M51 HDR images that I wanted to save to a folder called `m51_images_cleaned`, I would run the following command:
```
python3 image_review.py --outfolder m51_images_cleaned m51*HDR.fits
```

# Choosing Files
Although it is possible to type out each of your fits files, I highly recommend learning to
use wildcards. For example, SEO images currently (as of April 2022) look like this:
```
ngc2175_h-alpha_128.0s_bin1H_220406_042522_ordelore_seo_4_RAW.fits
```
Put another way, each fits file is named as follows:
```
{object name}_{filter name}_{exposure time}_{binning}{HDR gain}_{date}_{time}_{username}_{telescope}_{number}_{pipestep}.fits
```
If I wanted to select each RAW NGC2175 image, I could use `ngc2175*RAW.fits` as the filename. If I wanted
each FCAL h-alpha image of M51, I would use `m51_h-alpha*FCAL.fits`. If I wanted each file that was taken
on March 1st, 2022, I would use `*_220301_*.fits`. Each of these wildcards will search through the avilable
fits files in the current folder and return a list of them to image_filter.py. If you wanted to search through
other directories, you just need to append the directory name to the file: `/path/to/dir/*.fits`. You can also
use wildcards in directories to: `/path/to/dir/*/m51*.fits` will take each image of M51 in each of the directories.