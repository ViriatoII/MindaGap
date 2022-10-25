doc = ''' USAGE:   python mindagap.py  <PANORAMA.tif> <boxsize> <loopnum> --edges <True|False>

   Takes a single panorama image and fills the empty grid lines with neighbour-weighted values.
   Small box sizes yield limited results but work the best with a high loop number (like 20)
   
   Increase boxsize (s) to overcome bigger gaps 
   --edges is optional parameter to blur area around grid, for smoother transitions between tiles with different exposures (EXPERIMENTAL)
   
   
    27/06/2022
    Ricardo Guerreiro
    Resolve Biosciences'''


def fill_grids(img_array, box_size = 5, nloops = 1, edges = False):
    ''' Fills grid locations (pixel val 0) with values from neighbourhood through a gaussian kernel 
    Small box sizes yield limited results but work the best with a high loop number (like 20)'''

    # Grid coordinates and a copy of image
    grid_coords = img_array == 0
    im_copy = img_array.copy()

    # Make grid pixels have the minimum value of image (excluding 0 of grid)
    im_copy[grid_coords] = min(img_array[img_array > 0].flatten()) 
   
    if edges:
         # Kernel: here you should define how much the True "dilates"
         expansion = np.array([[True] *20] *20)

         # Expand/dilate grid   (https://python.tutorialink.com/how-to-expand-dilate-a-numpy-array/)
         # Convolution is not possible for bool values, so we convert to int and back. That works because bool(N) == True if N != 0.
         expanded_grid = convolve2d(grid_coords.astype(int), 
                                    expansion.astype(int), mode='same').astype(bool)
   
    # Create a blurred image and replace original grid positions by new blur value. Iterate
    for i in range(nloops):

        # Smooth edges as well (last loops only) if asked
        if edges and (i > nloops -5):
            blur_img = cv2.GaussianBlur(im_copy,(box_size,box_size), 0)  
            im_copy[expanded_grid] = blur_img[expanded_grid]     
            
        else: # Main condition
            blur_img = cv2.GaussianBlur(im_copy,(box_size,box_size), 0)   # Gaussian kernel
            im_copy[grid_coords] = blur_img[grid_coords]

    return(im_copy)


#################################################################################

                        
if __name__ == '__main__':
    import os,glob, argparse
    from scipy.signal import convolve2d
    import numpy as np
    import matplotlib.pyplot as plt

    #increase max allowed image size
    os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()

    import tifffile, cv2  
   

    parser = argparse.ArgumentParser(description="Takes a single panorama image and fills the empty grid lines with neighbour-weighted values" )
    parser.add_argument("input", help="Input tif/png file with grid lines to fill")
    parser.add_argument("s", nargs = '?', type=int, default = 5,  help="Box size for gaussian kernel (bigger better for big gaps but less accurate)")
    parser.add_argument("r", nargs = '?', type=int, default = 40, help="Number of rounds to apply gaussianBlur (more is better)")
    parser.add_argument("-e", '--edges', nargs = '?', default = False, help="Also smooth edges near grid lines")
    args=parser.parse_args()

    if args.s % 2 == 0:
        print("-s argument must be uneven number") ; exit()

    # Sanity checks of input
    pathname, extension = os.path.splitext(args.input)

    if os.path.exists(args.input) == False:
        print("Input file does not exist!") ; exit()

    # Read input as tif file or as png/
    if extension[1:4] == 'tif':
        img = tifffile.imread(args.input) 
    else:
        img = plt.imread(args.input)

    # Apply fill_grids function and write to file #####
        # Work on composite images or z-layered tiffs
    if len(img.shape) > 2: 
        layers = []
        for l in range(img.shape[0]):
            layers.append(fill_grids(img_array=img[l,:,:], box_size = args.s, nloops= args.r, edges = args.edges))
        img = np.array(layers )

    else: # Work on 2D images
        img = fill_grids(img_array=img, box_size = args.s, nloops= args.r, edges = args.edges)



    # Read input as tif file or as png/
    if extension[1:4] == 'tif':
        tifffile.imwrite(pathname + 'gridfilled' + extension, img)
    else:
        plt.imsave(pathname + 'gridfilled' + extension, img)
