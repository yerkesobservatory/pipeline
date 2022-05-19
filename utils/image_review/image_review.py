import sys
import argparse
import os
from turtle import back
from black import out
import numpy as np
import shutil
import pygame
import pygame.freetype
import matplotlib.colors as colors
import matplotlib.cm as cm
from astropy.stats import mad_std
import astropy.io.fits as fits
import yaml

# Constants
# If the image goes beyond the bounds of your screen, feel free to scale it down
SCREEN_SIZE = (800,800)

# If the filenames extned beyodn the bounds of your screen,
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
SCALING_FUNCS = [notebook_scaling, scaling_90, scaling_95, scaling_99, even_scaling]

def main(out_path, filenames, startindex, displaykeys, hduindex, displayhduname):
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
    
    while True:
        redraw = False
        fits_file = filenames[current_image]
        fits_hdulist = fits.open(fits_file)
        if hduindex < 0 or hduindex >= len(fits_hdulist):
            print("HDU index out of bounds")
            return
        fits_data = fits_hdulist[hduindex].data
        fits_header = fits_hdulist[hduindex].header
        
        # draw the fits image as the background
        fits_img = conv_to_cmap(fits_data, SCALING_FUNCS[scaling_index])
        display.fill('black')
        fits_buffer = pygame.image.frombuffer(fits_img.tobytes(), (fits_img.shape[1], fits_img.shape[0]), 'RGBA')
        display.blit(pygame.transform.scale(fits_buffer, SCREEN_SIZE), (0, 0))
        
        # filename and the image index
        filename_text = screen_font.render(f"{os.path.basename(filenames[current_image])} {current_image+1}/{len(filenames)}", bgcolor="white", fgcolor="black")[0]
        
        # fits keywords and other information
        info_row = ""
        if (displayhduname):
            info_row += f"{fits_header['EXTNAME']} | "
        if (os.path.exists(os.path.join(out_folder, os.path.basename(filenames[current_image])))):
            info_row += "Saved | "
        else:
            info_row += "Not Saved | "
        
        for key in displaykeys:
            if key in fits_header:
                info_row += f"{displaykeys[key]}: {fits_header[key]} | "
        file_info = screen_font.render(info_row, bgcolor="white", fgcolor="black")[0]
        
        # Report out background level & RMS (Add to keywords from SEP)
        
        # Display background subtracted image if it's an FCAL image
        
        display.blit(filename_text, (0,0))
        display.blit(file_info, (0,SCREEN_SIZE[1]-FONT_SIZE))
        
        
        # redraw the screen
        pygame.display.update()
        
        while not(redraw):
            for event in pygame.event.get():
                # if the user presses the quit button on the window (red button on Mac, [x] on Windows)
                # then quit the program
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()
                
                # event.key is only populated when the event is a key press
                if event.type == pygame.KEYDOWN:
                    # if the user presses escape, then quit
                    if event.key == pygame.K_ESCAPE:
                        pygame.quit()
                        sys.exit()
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
                        out_file = os.path.basename(fits_file)
                        shutil.copy(fits_file, os.path.join(out_path,out_file))
                    # if the user presses delete, then delete the image from the output folder
                    if event.key == pygame.K_BACKSPACE:
                        redraw = True
                        out_file = os.path.basename(fits_file)
                        # only delete the file if it already exists
                        if os.path.exists(os.path.join(out_path, out_file)):
                            os.remove(os.path.join(out_path, out_file))
                    # if the user presses a number, then change the colormap
                    if event.key in [pygame.K_1, pygame.K_2, pygame.K_3, pygame.K_4, pygame.K_5, pygame.K_6, pygame.K_7, pygame.K_8, pygame.K_9]:
                        new_scaling_idx = event.key - pygame.K_1
                        if new_scaling_idx < len(SCALING_FUNCS):
                            scaling_index = new_scaling_idx
                        redraw = True

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preview and filter FITS files')
    parser.add_argument("--outfolder", dest="out_folder", metavar="O", type=str, nargs=1, help="The folder used to save filtered images", required=False)
    parser.add_argument("--startindex", dest="startindex", metavar="O", type=int, nargs=1, help="The integer index of the start image (Default: 0)", required=False, default=0)
    parser.add_argument("--configfile", dest="config", metavar="O", type=str, nargs=1, help="The name of the configuration file (if it exists)", required=False)
    parser.add_argument("file", metavar="F", type=str, nargs=argparse.REMAINDER, help="The FITS files to be loaded and previewed")
    args = parser.parse_args()
    
    out_folder = None
    displaykeys = {}
    
    fits_files = args.file
    config_file = args.config
    startindex = 0
    hduindex = 0
    displayhduname = False
    # if the configuration file exists, then use it to form a guess for the flag values
    # if values are passed in as arguments, then they override the values in the config file
    if (config_file != None):
        config_stream = open(config_file[0], "r")
        yaml_file = yaml.load(config_stream, Loader=yaml.FullLoader)
        if 'out_folder' in yaml_file:
            out_folder = yaml_file['out_folder']
        if 'startindex' in yaml_file:
            startindex = yaml_file['startindex']
        if 'displaykeys' in yaml_file:
            displaykeys = yaml_file['displaykeys']
        if 'hduindex' in yaml_file:
            hduindex = yaml_file['hduindex']
        if 'displayhduname' in yaml_file:
            displayhduname = yaml_file['displayhduname']
            
    if args.out_folder != None:
        out_folder = args.out_folder[0]
    if type(args.startindex) == list:
        startindex = args.startindex[0]
    
    if out_folder == None:
        print("No output folder specified. Exiting...")
    elif fits_files == None or len(fits_files) == 0:
        print("No FITS files specified. Exiting...")
    else:
        if not(os.path.exists(out_folder)):
            print(f"{out_folder} is not a directory. Creating it now")
            os.mkdir(out_folder)
        
        main(out_folder, fits_files, startindex, displaykeys, hduindex, displayhduname)