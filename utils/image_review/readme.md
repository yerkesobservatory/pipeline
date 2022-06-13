# Installing
Just run `pip install -r requirements.txt` to install the necessary Python dependencies.

# Running
To run the image review script, run the following line in your terminal.
```
python3 image_review.py <input files>
```
If you want to use more options,
```
usage: image_review.py [-h] [--startindex O] [--config O] ...

Preview and filter FITS files

positional arguments:
  F               The FITS files to be loaded and previewed

optional arguments:
  -h, --help      show this help message and exit
  --startindex O  The integer index of the start image (Default: 0)
  --config O      The name of the configuration file
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
use wildcards in directories too: `/path/to/dir/*/m51*.fits` will take each image of M51 in each of the directories.

# Configuration Files and Saving Piperuns
To streamline the process of running this script, it is possible to update settings in a configuration file. This
avoids the need to type out long commands and speeds up the process of changing specific flags. Confiugrations files can be passed to the program using the `--config` flag.

Configuration files are stored as headers in configuration files. If you are using deltaconfigs for your personal versions of the SEO pipeline, I recommend appending this `image_filter` section to your deltaconfig file. Otherwise, you can use a text file that only contains the following section.

## Saving piperuns
Although it is possible to use this program to just review images, it is also possible to create a new piperun and use it to selectively process images. This can be useful if you want to stack images and want to avoid images that have lost tracking or images with poor seeing. In order to save a piperun, you MUST use a configuration file because the flags required to pass the piperun arguments would become unwieldy.

```
[image_filter]
    # If a file has multiple HDUs, this is the index of the HDU to load
    hduindex = 0

    # Decide whether to display the name of the HDU on screen
    displayhduname = True

    # Out of the loaded FITS files, start with the following index
    startindex = 0

    # Decide whether to save a new piperun that will allow you to further reduce images
    # by running `darepyperun.py piperun_filtered.txt`
    # These settings can be copied from a piperun you already use
    save_piperun = True
    pythonpath = "/Users/ordelore/Documents/SEO/content/pipeline/source"
    pipeconf = "/Users/ordelore/Documents/SEO/content/pipeline/config/pipeconf_SEO.txt /Users/ordelore/Documents/SEO/content/pipeline/config/lorenzo_config.txt"
    pipemode = "filter_second_half"
    logfile = "/Users/ordelore/Documents/SEO/content/pipelog.txt"
    loglevel = "DEBUG"
    # Make sure to change this line if you want to use a different output folder
    # if you are reducing a different target.
    outputfolder = "/Users/ordelore/Documents/SEO/content/data/m81"

    # to add more displayed fields, add them in this format
    # header name = "display text"
    # In the program, they will be displayed as
    # "display text: {value of header name}"
    [[displaykeys]]
        DEWTEM1 = "Sensor Temp"
        RHALF = "Half radius"
```