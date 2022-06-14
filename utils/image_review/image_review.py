'''
Title: image_review.py
Description: This program will allow the user to review a batch set of FITS
    images and allows them to potentially output those filenames into a piperun
    that can be interpreted by the SEO pipeline for further processing
Author: Lorenzo Orders
'''

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
import webbrowser
from astropy.visualization import ZScaleInterval

# Constants
# If the image goes beyond the bounds of your screen, feel free to scale it down
SCREEN_SIZE = (800, 800)

# If the filenames extend beyond the bounds of your screen, you can scale down
# the size of the text
FONT_SIZE = 18

# When multiple lines are rendered onscreen, at once, this determines the spacing
# between each line.
TEXT_PADDING = 2

# these are the different colormaps the program can display. If you want to add more,
# visit https://matplotlib.org/stable/gallery/color/colormap_reference.html
COLORMAP_NAMES = ['gist_yarg', 'gist_gray', 'turbo', 'gist_ncar', 'Set1']


def conv_to_cmap(image, scaling_func, cmapidx):
    """
    conv_to_cmap will take a 2D array, a scaling function, and an index for the
    colormap to use it will return an RGB image where the scaling function and
    colormap are used to convert each value to an RGB color

    Parameters
    ----------
    image : ndarray
        The 2D array of scalar values to convert to an RGB image
    scaling_func : function (image : ndarray) -> (vmin : float, vmax : float)
        The function that will calculate vmin and vmax for the image
    cmapidx : int
        The index of the matplotlib colormap to use

    Returns
    -------
    rgb_image : ndarray
        A 3D array of RGB values. The first two dimensions are the same as the
        input image.
    """
    vmin, vmax = scaling_func(image)
    cmap = cm.get_cmap(COLORMAP_NAMES[cmapidx])
    norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    return mapper.to_rgba(image, bytes=True)


def notebook_scaling(image):
    """
    notebook scaling will take a 2D ndarray and return the vmin and vmax that
    should form the minimum and maximum color values. Any values beyond those
    bounds will be clamped to the nearest bound. Each of these scaling functions
    have the same function signature, so only this function will have a docstring.
    If you create your own function, it is vital to match the function signature
    exactly.
    
    Parameters
    ----------
    image : ndarray
        The 2D array of scalar values to form the basis of the vmin/vmax calculation

    Returns
    -------
    vmin : float
        The minimum value that should be used to form the color map
    vmax : float
        The maximum value that should be used to form the color map
    """
    med = np.nanmedian(image)
    mad = mad_std(image, ignore_nan=True)
    vmin = med - mad
    vmax = med + 3 * mad
    return vmin, vmax


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
    vmin = med - 3 * mad
    vmax = med + 3 * mad
    return vmin, vmax


def zscale_scaling(image):
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(image)
    return vmin, vmax


# if you make your own scaling function, remember to add it to the list
SCALING_FUNCS = [notebook_scaling, scaling_90,
                 scaling_95, scaling_99, even_scaling, zscale_scaling]


def render_multiline_text(text, font_renderer):
    """
    render_multiline_text will take a string and a font renderer and return
    a new pygame surface where the string is displayed and none of the text
    goes off the screen
    
    Parameters
    ----------
    text : str
        The string to display
    font_renderer : pygame.freetype.Font
        The font renderer to use to render the text
    
    Returns
    -------
    text_surface : pygame.Surface
        The surface containing the text
    text_render_pos : tuple
        Assuming that the text is rendered at the bottom of the screen, this
        is the position to render text_surface at
    """
    words = text.split(' | ')
    line_len = 0
    screen_width = SCREEN_SIZE[0]
    screen_height = SCREEN_SIZE[1]
    # surface to render the text onto
    text_surface = pygame.Surface(SCREEN_SIZE)
    text_surface.fill((255, 255, 255))
    # the position to render the text surface on the main screen
    _, size_info = font_renderer.render(text[0])
    text_height = size_info.height

    text_render_pos = [0, screen_height - text_height - TEXT_PADDING]
    lines = 0
    for word in words:
        if len(word) == 0:
            continue
        word_surface, word_rect = font_renderer.render(
            word + " | ", bgcolor="white", fgcolor="black")
        line_len += word_rect.width
        # if rendering this word would have exceeded the length of the
        # screen, then go to a new line
        if line_len > screen_width:
            text_render_pos[1] -= word_rect.height
            line_len = word_rect.width
            lines += 1
        text_pen = (line_len - word_rect.width,
                    (text_height + TEXT_PADDING) * lines)
        text_surface.blit(word_surface, text_pen)
    return text_surface, text_render_pos


def main(filenames, displaykeys, hduindex, displayhduname):
    pygame.init()
    display = pygame.display.set_mode(SCREEN_SIZE)
    pygame.display.set_caption("Image Review")

    current_image = 0
    pygame.freetype.init()

    screen_font = pygame.freetype.SysFont("Arial", size=FONT_SIZE)

    scaling_index = 0
    colormap_index = 0
    saving_files = set()

    while True:
        redraw = False
        fits_file = filenames[current_image]
        fits_hdulist = fits.open(fits_file)

        # load a specific hdu
        if hduindex < 0 or hduindex >= len(fits_hdulist):
            print("HDU index out of bounds")
            return saving_files

        fits_data = fits_hdulist[hduindex].data
        fits_header = fits_hdulist[hduindex].header
        primary_header = fits_hdulist[0].header
        fits_hdulist.close()

        # check for a binary table
        if 'XTENSION' in fits_header and fits_header['XTENSION'] == 'BINTABLE':
            print("Binary table detected. Please use a FITS image")
            return saving_files

        # draw the fits image
        fits_img = conv_to_cmap(
            fits_data, SCALING_FUNCS[scaling_index], colormap_index)
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
            # if there are two values from the config file, then one of them defines
            # the number of decimal places to display
            if type(displaykeys[key]) == type([]):
                displaytext = displaykeys[key][0]
                display_decimals = int(displaykeys[key][1])
                if key in fits_header:
                    info_row += f"{displaytext}: {round(fits_header[key], display_decimals)} | "
                elif key in primary_header:
                    info_row += f"{displaytext}: {round(primary_header[key], display_decimals)} | "
            else:
                if key in fits_header:
                    info_row += f"{displaykeys[key]}: {fits_header[key]} | "
                elif key in primary_header:
                    info_row += f"{displaykeys[key]}: {primary_header[key]} | "
        file_info_text, file_info_pos = render_multiline_text(
            info_row, screen_font)

        # draw information to the screen
        display.blit(filename_text, (0, 0))
        display.blit(file_info_text, file_info_pos)

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
                    # if the user presses space, and the image is not currently set to be saved,
                    # then save the image. Otherwise, remove the image from the list of images to process
                    if event.key == pygame.K_SPACE:
                        redraw = True
                        if filenames[current_image] in saving_files:
                            saving_files.remove(filenames[current_image])
                        else:
                            saving_files.add(filenames[current_image])
                    # if the user presses a number, then change the colormap
                    if event.key in [pygame.K_1, pygame.K_2, pygame.K_3, pygame.K_4, pygame.K_5, pygame.K_6, pygame.K_7, pygame.K_8, pygame.K_9]:
                        new_scaling_idx = event.key - pygame.K_1
                        if new_scaling_idx < len(SCALING_FUNCS):
                            scaling_index = new_scaling_idx
                        redraw = True
                    # if the user presses 'h', then open readme.md in a web browser
                    if event.key == pygame.K_h:
                        webbrowser.open(
                            'https://github.com/yerkesobservatory/pipeline/blob/dev/utils/image_review/readme.md')

                    # advance the colormap index
                    if event.key == pygame.K_COMMA and colormap_index > 0:
                        colormap_index -= 1
                        redraw = True
                    # retreat the colormap index
                    if event.key == pygame.K_PERIOD and colormap_index < len(COLORMAP_NAMES) - 1:
                        colormap_index += 1
                        redraw = True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Preview and filter FITS files')
    parser.add_argument("--config", dest="config", metavar="O", type=str, nargs=1,
                        help="The name of the configuration file", required=False)
    parser.add_argument("file", metavar="F", type=str, nargs=argparse.REMAINDER,
                        help="The FITS files to be loaded and previewed")
    args = parser.parse_args()

    displaykeys = {}

    fits_files = args.file
    config_file = args.config
    hduindex = 0
    displayhduname = False
    config_obj = None
    saving_piperun = False

    # if the configuration file exists, then use it to form a guess for the flag values
    # if values are passed in as arguments, then they override the values in the config file
    if (config_file != None):
        config_obj = configobj.ConfigObj(config_file[0])['image_filter']
        if 'displaykeys' in config_obj:
            displaykeys = config_obj['displaykeys']
        if 'hduindex' in config_obj:
            hduindex = int(config_obj['hduindex'])
        if 'displayhduname' in config_obj:
            displayhduname = config_obj['displayhduname'] == "True"

    saving_piperun = config_obj and 'save_piperun' in config_obj and config_obj[
        'save_piperun'] == "True"
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
        fits_files, displaykeys, hduindex, displayhduname)

    # saves a piperun so the user can run 'darepyperun.py piperun_filtered' to continue processing their images
    if saving_piperun and len(files_to_process) > 0:
        if 'piperun_name' in config_obj:
            outfilename = config_obj['piperun_name']
        else:
            outfilename = "piperun_filtered.txt"
        piperun_out = f"pythonpath = {config_obj['pythonpath']}\n"
        config_files = config_obj['pipeconf'].replace(' ', '\n')
        piperun_out += f"pipeconf = {config_files}\n"
        piperun_out += f"pipemode = {config_obj['pipemode']}\n"
        piperun_out += f"logfile = {config_obj['logfile']}\n"
        piperun_out += f"loglevel = {config_obj['loglevel']}\n"
        piperun_out += f"outputfolder = {config_obj['outputfolder']}\n"
        files_out = '\n'.join(files_to_process)
        piperun_out += f"inputfiles = {files_out}\n"
        with open(outfilename, 'w') as f:
            f.write(piperun_out)
    sys.exit()
