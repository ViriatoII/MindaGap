doc = ''' USAGE:   python mindagap.py  <PANORAMA.tif> <boxsize> <loopnum> -xt <Xtilesize> -yt <Ytilesize> --edges <True|False> 

   Takes a single panorama image and fills the empty grid lines with neighbour-weighted values.
   
   Small boxsize yields limited results but works the best with a high loop number (like 20)
   Increase boxsize to overcome bigger gaps 
   <tilesize> is optional parameter to deal with adjacent tiles without gaps but where there are visible different exposures  (EXPERIMENTAL) 
   
   --edges is optional parameter to blur area around grid, for smoother transitions between tiles with different exposures (EXPERIMENTAL)
   
   
    27/06/2022
    Ricardo Guerreiro
    Resolve Biosciences'''      

    
def read_img(path_file):
    ''' Reads a tiff/png image into a numpy array'''
    
    pathname, extension = os.path.splitext(path_file)

    try:
        if not os.path.exists(path_file):
            print(f"Input file does not exist: \n {path_file}") 
            
            return(0)

        # Read input as tif file or as png/jpg
        if extension[1:4] == 'tif':
            img = tifffile.imread(path_file) 
        else:
            img = cv2.imread(path_file,cv2.IMREAD_UNCHANGED) 

        return(img)
    except: 
        print(f"Input file invalid: \n {path_file}") 


def fill_grids(img_array, box_size = 5, nloops = 1, edges = 0, Xtilesize = None):
    ''' Fills grid locations (pixel val 0) with values from neighbourhood through a gaussian kernel 
    Small box sizes yield limited results but work the best with a high loop number (like 20)'''

    # Grid coordinates and a copy of image
    grid_coords = img_array == 0
    im_copy = img_array.copy()   
    edges = int(edges)
   
    # Make grid pixels have the minimum value of image (excluding 0 of grid)
    im_copy[grid_coords] = min(img_array[img_array > 0].flatten()) 

    if Xtilesize: 
       # First iteration going through vertical lines 
       for yjump in range(Ytilesize, panoYmax, Ytilesize):
           ymin,ymax = yjump -1 ,yjump + 2
           grid_coords [ymin:ymax,:] = True 

       for xjump in range(Xtilesize, panoXmax, Xtilesize):
           xmin,xmax = xjump -1 ,xjump + 2
           grid_coords [:,xmin:xmax] = True 
    print("Test 2")
   

    if edges > 0:
         # Kernel: here you should define how much the True "dilates"

         expansion = np.array([[True] *edges] *edges)

         # Expand/dilate grid   (https://python.tutorialink.com/how-to-expand-dilate-a-numpy-array/)
         # Convolution is not possible for bool values, so we convert to int and back. That works because bool(N) == True if N != 0.
         expanded_grid = convolve2d(grid_coords.astype(int), 
                                    expansion.astype(int), mode='same').astype(bool)
   
    # Create a blurred image and replace original grid positions by new blur value. Iterate
    for i in track(range(nloops),description='[green]Applying Gaussian Blur to gridlines ...'):

        # Smooth edges as well (last loops only) if asked
        if edges and (i > nloops -5):
            blur_img =  cv2.medianBlur(im_copy,box_size,cv2.BORDER_DEFAULT) 
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
    from rich.progress import track

    #increase max allowed image size
    os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
    
    version_number = "0.0.4"

    import tifffile,  cv2  

    parser = argparse.ArgumentParser(description="Takes a single panorama image and fills the empty grid lines with neighbour-weighted values" )
    parser.add_argument("input", help="Input tif/png file with grid lines to fill")
    parser.add_argument("-s", "--sizekernel", nargs = '?', type=int, default = 5,  help="Box size for gaussian kernel (bigger better for big gaps but less accurate)")
    parser.add_argument("-r", "--rounds",    nargs = '?', type=int, default = 40, help="Number of rounds to apply gaussianBlur (more is better)")
    parser.add_argument("-xt", "--Xtilesize", nargs = '?', type=int, default = None,  help="Tile size (distance between gridlines) on X axis")
    parser.add_argument("-yt", "--Ytilesize", nargs = '?', type=int, default = None,  help="Tile size (distance between gridlines) on Y axis")

    parser.add_argument("-e", '--edges', nargs = '?', default = 0, help="Also smooth edges near grid lines")
    parser.add_argument("-v", '--version', action='store_true', default = False, help="Print version number.")
    args=parser.parse_args()


    if args.sizekernel % 2 == 0:
        print("-s argument must be uneven number") ; exit()
        
    if args.version:        print(version_number) ; exit()

    # Sanity checks of input
    pathname, extension = os.path.splitext(args.input)

    if os.path.exists(args.input) == False:
        print("Input file does not exist!") ; exit()

    # Inputs
    Xtilesize = args.Xtilesize #2144
    Ytilesize = args.Xtilesize if args.Ytilesize == None else args.Ytilesize 


    # Read input as tif file or as png/jpg
    img = read_img(args.input) 
    panoYmax,panoXmax  = img.shape [-2:]
    
    print("Test 1")

    # Apply fill_grids function and write to file #####
        # Work on composite images or z-layered tiffs
    if len(img.shape) > 2: 
        layers = []
        for l in range(img.shape[0]):

            layers.append(fill_grids(img_array=img[l,:,:], box_size = args.sizekernel, nloops= args.rounds, edges = args.edges))
        img = np.array(layers )

    else: # Work on 2D images

        img = fill_grids(img_array=img, box_size = args.sizekernel, nloops= args.rounds, edges = args.edges)



    # Save as tif file or as png/

    #cv2.imwrite(pathname + '_gridfilled' + extension, img)

    if 'tif' in extension:
        tifffile.imwrite(pathname + '_gridfilled' + extension, img)
    else:
        cv2.imwrite(pathname + '_gridfilled' + extension, img)
