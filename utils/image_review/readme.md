# Installing
The easiest way to download this program is to download it as part of the
pipeline. If you are not using the pipeline, or find this tool useful for
other purposes, then you can download it by clicking [this link](https://raw.githubusercontent.com/yerkesobservatory/pipeline/dev/utils/image_review/image_review.py) and
then save the file to a folder on your computer. Make sure to also download
the [requirements.txt](https://raw.githubusercontent.com/yerkesobservatory/pipeline/dev/utils/image_review/requirements.txt) file as well and save it to the same location.

## Dependencies
To install the necessary Python packages, `cd` to the folder where image_review
is located and run `pip install -r requirements.txt` to install the necessary
Python dependencies.

# Running
To run the image review script, run the following line in your terminal.
```
python3 image_review.py <input files>
```
If you want to use more options,
```
usage: image_review.py [-h] [--config O] ...

Preview and filter FITS files

positional arguments:
  F           The FITS files to be loaded and previewed

optional arguments:
  -h, --help  show this help message and exit
  --config O  The name of the configuration file
  ```

# Keybindings
<kbd>Esc</kbd>: Close the program and save a piperun if necessary

<kbd>&#8592;</kbd>-<kbd>&#8594;</kbd>: Step through each image

<kbd>0</kbd>-<kbd>9</kbd>: Switch between different color scales

<kbd>Space</kbd>: Toggle between saving and not saving the image

<kbd>h</kbd>: Opens this readme in a web browser

<kbd>,</kbd>: Go back through the list of colormaps

<kbd>.</kbd>: Go forward through the list of colormaps

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
If I wanted to select each RAW NGC2175 image, I could use `ngc2175*RAW.fits` as
the filename. If I wanted each FCAL h-alpha image of M51, I would use
`m51_h-alpha*FCAL.fits`. If I wanted each file that was taken on March 1st, 2022,
I would use `*_220301_*.fits`. Each of these wildcards will search through the
available fits files in the current folder and return a list of them to image_filter.py.
If you wanted to search through other directories, you just need to append the directory
name to the file: `/path/to/dir/*.fits`. You can also use wildcards in directories,
too: `/path/to/dir/*/m51*.fits` will take each image of M51 in each of the directories.

# Configuration Files and Saving Piperuns
To streamline the process of running this script, it is possible to update
settings in a configuration file. This avoids the need to type out long commands
and speeds up the process of changing specific flags. Configurations files can be
passed to the program using the `--config` flag.

Configuration settings for image_review are stored as headers in configuration files. If you are using deltaconfigs in your personal versions of the SEO pipeline, I recommend appending a `image_filter` section to your deltaconfig file. Otherwise, you can use a text file that only contains an [image_filter] section.

## Saving piperuns
Although it is possible to use this program to just review images, it is also possible to create a new piperun and use it to selectively process images. This can be useful if you want to stack images and want to avoid images that have lost tracking or images with poor seeing. In order to save a piperun, you MUST use a configuration file because the flags required to pass the piperun arguments would become unwieldy.

```
[image_filter]
    # If a file has multiple HDUs, this is the index to load
    hduindex = 1

    # Decide whether to display the name of the HDU
    displayhduname = True

    # the name of the file to save the piperun to
    piperun_name = "piperun_filtered.txt"

    # if save_piperun is True, and all the below variables are defined, then this program will output a new piperun
    save_piperun = True
    pythonpath = "/Users/ordelore/Documents/SEO/content/pipeline/source"
    pipeconf = "/Users/ordelore/Documents/SEO/content/pipeline/config/pipeconf_SEO.txt /Users/ordelore/Documents/SEO/content/pipeline/config/lorenzo_config.txt"
    pipemode = "filter_second_half"
    logfile = "/Users/ordelore/Documents/SEO/content/pipelog.txt"
    loglevel = "DEBUG"
    outputfolder = "/Users/ordelore/Documents/SEO/content/data/m3"

    # to add more displayed fields, add them in this format
    # header name = "display text", decimal places
    # In the program, they will be displayed as "display text: {value of header name truncated to the number of decimal places}"
    [[displaykeys]]
        DEWTEM1 = "Sensor Temp"
        RHALF = "Half radius", 3
        RHALFSTD = "Half radius std", 3
```