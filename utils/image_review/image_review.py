# Code written by Lorenzo Orders in Winter and Spring Quarters of 2022.

import sys
import argparse
import os
import numpy as np
import pygame
import pygame.freetype
import matplotlib.colors as colors
import matplotlib.cm as cm
from astropy.stats import mad_std
import astropy.io.fits as fits
import configobj

# Constants
# If the image goes beyond the bounds of your screen, feel free to scale it down
SCREEN_SIZE = (800, 800)

# If the filenames extned beyond the bounds of your screen,
# you can scale down the size of the text
FONT_SIZE = 18

# This function converts a 2D array of scalar values to a 2D array of RGB values
# min/max is borrowed from Al's imgplot function
# scaling_func MUST be a function that takes in the image data and returns the
# minimum and maximum values for the colormap


def conv_to_cmap(image, scaling_func):
    # vmin = med-mad
    # vmax = med + 3*mad
    vmin, vmax = scaling_func(image)
    cmap = cm.get_cmap('binary')
    norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    return mapper.to_rgba(image, bytes=True)

# based on Al's notebook scaling


def notebook_scaling(image):
    med = np.nanmedian(image)
    mad = mad_std(image, ignore_nan=True)
    return med-mad, med + 3*mad


def scaling_90(image):
    vmax = np.quantile(image.flatten(), 0.9)
    vmin = np.quantile(image.flatten(), 0.1)
    return vmin, vmax


def scaling_95(image):
    vmax = np.quantile(image.flatten(), 0.95)
    vmin = np.quantile(image.flatten(), 0.05)
    return vmin, vmax


def scaling_99(image):
    vmax = np.quantile(image.flatten(), 0.99)
    vmin = np.quantile(image.flatten(), 0.01)
    return vmin, vmax


def even_scaling(image):
    med = np.nanmedian(image)
    mad = mad_std(image, ignore_nan=True)
    return med-3*mad, med + 3*mad


# if you make your own scaling function, remember to add it to the list
SCALING_FUNCS = [notebook_scaling, scaling_90,
                 scaling_95, scaling_99, even_scaling]


def main(filenames, startindex, displaykeys, hduindex, displayhduname):
    pygame.init()
    display = pygame.display.set_mode(SCREEN_SIZE)
    pygame.display.set_caption("Image Review")
    if startindex < 0:
        startindex = 0
    if startindex >= len(filenames):
        startindex = len(filenames)-1

    current_image = startindex
    pygame.freetype.init()

    screen_font = pygame.freetype.SysFont("Arial", size=FONT_SIZE)

    scaling_index = 0
    saving_files = set()

    while True:
        redraw = False
        fits_file = filenames[current_image]
        fits_hdulist = fits.open(fits_file)

        # load a specific hdu
        if hduindex < 0 or hduindex >= len(fits_hdulist):
            print("HDU index out of bounds")
            return
        fits_data = fits_hdulist[hduindex].data
        fits_header = fits_hdulist[hduindex].header
        fits_hdulist.close()

        # draw the fits image as the background
        fits_img = conv_to_cmap(fits_data, SCALING_FUNCS[scaling_index])
        display.fill('black')
        fits_buffer = pygame.image.frombuffer(
            fits_img.tobytes(), (fits_img.shape[1], fits_img.shape[0]), 'RGBA')
        display.blit(pygame.transform.scale(fits_buffer, SCREEN_SIZE), (0, 0))

        # filename and the image index
        filename_text = screen_font.render(
            f"{os.path.basename(filenames[current_image])} {current_image+1}/{len(filenames)}", bgcolor="white", fgcolor="black")[0]

        # fits keywords and other information
        info_row = ""
        if (displayhduname):
            info_row += f"{fits_header['EXTNAME']} | "
        if (filenames[current_image] in saving_files):
            info_row += "Saved | "
        else:
            info_row += "Not Saved | "

        for key in displaykeys:
            if key in fits_header:
                info_row += f"{displaykeys[key]}: {fits_header[key]} | "
        file_info = screen_font.render(
            info_row, bgcolor="white", fgcolor="black")[0]

        # draw information to the screen
        display.blit(filename_text, (0, 0))
        display.blit(file_info, (0, SCREEN_SIZE[1]-FONT_SIZE))

        # redraw the screen
        pygame.display.update()

        while not(redraw):
            for event in pygame.event.get():
                # if the user presses the quit button on the window (red button on Mac, [x] on Windows)
                # then quit the program
                if event.type == pygame.QUIT:
                    pygame.quit()
                    return saving_files

                # event.key is only populated when the event is a key press
                if event.type == pygame.KEYDOWN:
                    # if the user presses escape, then quit
                    if event.key == pygame.K_ESCAPE:
                        pygame.quit()
                        return saving_files
                    # if the user presses left, then go to the previous fits file
                    if event.key == pygame.K_LEFT and current_image > 0:
                        current_image -= 1
                        redraw = True
                    # if the user presses right, then go to the next fits file
                    if event.key == pygame.K_RIGHT and current_image < len(filenames)-1:
                        current_image += 1
                        redraw = True
                    # if the user presses space, then save the image to the output folder
                    if event.key == pygame.K_SPACE:
                        redraw = True
                        saving_files.add(filenames[current_image])
                    # if the user presses delete, then delete the image from the output folder
                    if event.key == pygame.K_BACKSPACE:
                        redraw = True
                        # only delete the file if it already exists
                        if (filenames[current_image] in saving_files):
                            saving_files.remove(filenames[current_image])
                    # if the user presses a number, then change the colormap
                    if event.key in [pygame.K_1, pygame.K_2, pygame.K_3, pygame.K_4, pygame.K_5, pygame.K_6, pygame.K_7, pygame.K_8, pygame.K_9]:
                        new_scaling_idx = event.key - pygame.K_1
                        if new_scaling_idx < len(SCALING_FUNCS):
                            scaling_index = new_scaling_idx
                        redraw = True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Preview and filter FITS files')
    parser.add_argument("--startindex", dest="startindex", metavar="O", type=int, nargs=1,
                        help="The integer index of the start image (Default: 0)", required=False, default=0)
    parser.add_argument("--config", dest="config", metavar="O", type=str, nargs=1,
                        help="The name of the configuration file", required=False)
    parser.add_argument("file", metavar="F", type=str, nargs=argparse.REMAINDER,
                        help="The FITS files to be loaded and previewed")
    args = parser.parse_args()

    displaykeys = {}

    fits_files = args.file
    config_file = args.config
    startindex = 0
    hduindex = 0
    displayhduname = False
    config_obj = None
    
    # if the configuration file exists, then use it to form a guess for the flag values
    # if values are passed in as arguments, then they override the values in the config file
    if (config_file != None):
        config_obj = configobj.ConfigObj(config_file[0])['image_filter']
        print(config_obj)
        if 'startindex' in config_obj:
            startindex = int(config_obj['startindex'])
        if 'displaykeys' in config_obj:
            displaykeys = config_obj['displaykeys']
        if 'hduindex' in config_obj:
            hduindex = int(config_obj['hduindex'])
        if 'displayhduname' in config_obj:
            displayhduname = bool(config_obj['displayhduname'])
            
    if type(args.startindex) == list:
        startindex = args.startindex[0]

    saving_piperun = config_obj and 'save_piperun' in config_obj and bool(
        config_obj['save_piperun'])

    # if the user wants to save a piperun and the necessary flags are not present, then
    # do not even let them process images. This saves the time the user would have spent processing
    # the images only to find out they missed a single flag.
    if (saving_piperun and not('pythonpath' in config_obj and 'pipeconf' in config_obj and
                                'pipemode' in config_obj and 'logfile' in config_obj and
                                'loglevel' in config_obj and 'outputfolder' in config_obj)):
        print("Missing flags to save a piperun!")
        print("Please ensure that your configuration file contains \'pythonpath\', \'pipeconf\', \'pipemode\', \'logfile\', \'loglevel\', and \'outputfolder\' entries")
        sys.exit()
    
    elif fits_files == None or len(fits_files) == 0:
        print("No FITS files specified. Exiting...")
        parser.print_help()
        sys.exit()

    files_to_process = []
    files_to_process = main(
        fits_files, startindex, displaykeys, hduindex, displayhduname)

    # saves a piperun so the user can run 'darepyperun.py piperun_filtered' to continue processing their images
    if 'save_piperun' in config_obj and bool(config_obj['save_piperun']):
        piperun_out = f"pythonpath = {config_obj['pythonpath']}\n"
        config_files = config_obj['pipeconf'].replace(' ', '\n')
        piperun_out += f"pipeconf = {config_files}\n"
        piperun_out += f"pipemode = {config_obj['pipemode']}\n"
        piperun_out += f"logfile = {config_obj['logfile']}\n"
        piperun_out += f"loglevel = {config_obj['loglevel']}\n"
        piperun_out += f"outputfolder = {config_obj['outputfolder']}\n"
        files_out = '\n'.join(files_to_process)
        piperun_out += f"inputfiles = {files_out}\n"
        with open(f"piperun_filtered.txt", 'w') as f:
            f.write(piperun_out)
            print(piperun_out)
    sys.exit()
