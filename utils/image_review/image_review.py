import sys
import argparse
import os
from tracemalloc import start
import numpy as np
import shutil
import pygame
import pygame.freetype
import matplotlib.colors as colors
import matplotlib.cm as cm
from astropy.stats import mad_std
import astropy.io.fits as fits


# Constants
# If the image goes beyond the bounds of your screen, feel free to scale it down
SCREEN_SIZE = (800,800)

# If the filenames extned beyodn the bounds of your screen,
# you can scale down the size of the text
FONT_SIZE = 18

# This function converts a 2D array of scalar values to a 2D array of RGB values
# min/max is borrowed from Al's imgplot function
def conv_to_cmap(image):
    med = np.nanmedian(image)
    mad = mad_std(image, ignore_nan=True)
    vmin = med-mad
    vmax = med + 3*mad
    cmap = cm.get_cmap('binary')
    norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    return mapper.to_rgba(image, bytes=True)
    

def main(out_path, filenames, startindex):
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
    saved_message = screen_font.render("Saved", fgcolor="black", bgcolor="white")[0]
    deleted_message = screen_font.render("Not Saved", fgcolor="black", bgcolor="white")[0]
    while True:
        redraw = False
        fits_file = filenames[current_image]
        fits_data = fits.getdata(fits_file)
        fits_header = fits.getheader(fits_file)
        
        # draw the fits image as the background
        fits_img = conv_to_cmap(fits_data)
        display.fill('black')
        fits_buffer = pygame.image.frombuffer(fits_img.tobytes(), (fits_img.shape[1], fits_img.shape[0]), 'RGBA')
        display.blit(pygame.transform.scale(fits_buffer, SCREEN_SIZE), (0, 0))
        
        # filename and the image index
        filename_text = screen_font.render(f"{os.path.basename(filenames[current_image])} {current_image+1}/{len(filenames)}", bgcolor="white", fgcolor="black")[0]
        
        # image temperature
        info_row = ""
        if 'DEWTEM1' in fits_header:
            info_row += f"Sensor temp: {fits_header['DEWTEM1']}°C"
        if 'RHALF' in fits_header:
            info_row += f"  Half radius: {fits_header['RHALF']}px"
        file_temp = screen_font.render(info_row, bgcolor="white", fgcolor="black")[0]
        
        # Report out background level & RMS (Add to keywords from SEP)
        
        # Display background subtracted image if it's an FCAL image
        
        display.blit(filename_text, (0,0))
        display.blit(file_temp, (0,SCREEN_SIZE[1]-FONT_SIZE))
        
        # draw the saved message if the image exists in the output directory
        if (os.path.exists(os.path.join(out_folder, os.path.basename(filenames[current_image])))):
            display.blit(saved_message, (0,SCREEN_SIZE[1]-FONT_SIZE*2))
        else:
            display.blit(deleted_message, (0,SCREEN_SIZE[1]-FONT_SIZE*2))
        
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

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preview and filter FITS files')
    parser.add_argument("--outfolder", dest="out_folder", metavar="O", type=str, nargs=1, help="The folder used to save filtered images")
    parser.add_argument("--startindex", dest="startindex", metavar="O", type=int, nargs=1, help="The integer index of the start image (Default: 0)", required=False, default=0)
    parser.add_argument("file", metavar="F", type=str, nargs=argparse.REMAINDER, help="The FITS files to be loaded and previewed")
    args = parser.parse_args()
    out_folder = args.out_folder[0]
    fits_files = args.file
    if type(args.startindex) == list:
        startindex = args.startindex[0]
    else:
        startindex = args.startindex
    if not(os.path.exists(out_folder)):
        print(f"{out_folder} is not a directory. Creating it now")
        os.mkdir(out_folder)
    main(out_folder, fits_files, startindex)