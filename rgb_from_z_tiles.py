doc = ''' 
    USAGE:   python rgb_from_z_tiles.py  -b <DAPI.tiff> -r <WGA.tiff> -g <constructive.tiff>   -z 2,4-6,9,14:17,20

   Takes 3 input multi-stack tif files and creates RGB composite png file for each requested z-layer.
   If red and/or green layer are not provided, empty channels are outputed

  Mandatory named arguments:
  -b BLUE, --blue BLUE  DAPI tif
  -z [Z], --z [Z]       The desired Z layers to keep. Accepts individual layer numbers, slices or both, separated by commas (no spaces)

  Optional arguments:
  -outdir [OUT], --out [OUT]
                        Output directory to save RGB tif files
  -r [RED], --red [RED]
                        WGA (cellwall stain) tif
  -g [GREEN], --green [GREEN]
                        Constructive signal tif

  -gp GREEN_PADDING, --green_padding GREEN_PADDING
                        Add X green layers around asked Z

  -bp BLUE_PADDING, --blue_padding BLUE_PADDING
                        Average X blue layers around asked Z layer

  -corr CORRECT_ILUM, --correct_ilum CORRECT_ILUM
                        Correct ilumination for Red channel. 1.0 to 3.0 for fine tunning. 0 to turn off.



TODO: Clip maximum WGA and DAPI? -> Decreases sensitivity in overexposed areas but might increase sensitivity or less exposed
        - Change threshold system?  ---->  pcent = np.percentile(r,1) ;  rblur [r < pcent + 50 ] *= 2 

    13/07/2022
    Ricardo Guerreiro
    Resolve Biosciences'''


#################################################################################

                        
if __name__ == '__main__':
    import os,glob, argparse
    import numpy as np
    import tifffile as tif
    import cv2 

    # Parse input arguments
    parser = argparse.ArgumentParser(description="Reads 3 3D tif files, extracts desired z layer and creates 3-channel RGB tiff image" )
    #parser.add_argument("indir", help="Input directory with 3D tif files")
    parser.add_argument("-outdir", "--out", nargs = '?',  help= "Output directory to save RGB tif files")
    parser.add_argument('-r', '--red',   nargs = '?', default = None,  help="WGA (cellwall stain) tif") 
    parser.add_argument('-g', '--green', nargs = '?', default = None,  help="Constructive signal tif")  
    parser.add_argument('-gp', '--green_padding', type=int, default = 4, help="Add X green layers around asked Z layer")  
    parser.add_argument('-bp', '--blue_padding',  type=int, default = 2, help="Average X blue layers around asked Z layer")  
    parser.add_argument('-corr', '--correct_ilum', type=float, default = 2, help="Correct ilumination for Red channel. 1.0 to 3.0 for fine tunning. 0 to turn off.")  
    requiredNamed = parser.add_argument_group('Mandatory named arguments')
    requiredNamed.add_argument('-b', '--blue',  help="DAPI tif", required = True) 
    requiredNamed.add_argument('-z', '--z', nargs = '?',  help="The desired Z layers to keep")   # type=int, default = 5,
    args=parser.parse_args()


    outname = os.path.basename(args.blue).split('_')[0] 

    # Make sure outdir is absolute path
    if args.out:
        args.out = os.path.abspath(args.out)
        outname = os.path.join(args.out, outname)

    # Read mandatory Blue argument
    blue  = tif.imread(args.blue )        # '../../../13_rois_for_demo/32770-Rice-Benfey-slide1_W6A2_P1X_R8_RO-Channel3_04032022-10-25-26.ome.tiff')
    xend,yend = blue.shape[-2:]

    # Read greed and red tif files or create empty array if they are not provided
    if args.green:
        green = tif.imread(args.green)       #'../../../13_rois_for_demo/bc_32770-Rice-Benfey-slide1_W6A2_P1X_constructive.tif')    
    else:
        green = np.zeros(blue.shape)     

    if args.red:
        red   = tif.imread(args.red )         #'../../../13_rois_for_demo/32770-Rice-Benfey-slide1_W6A2_P1X_R9_*.ome.tiff') 
    
    else:
        red = np.zeros(blue.shape)   


    if args.z:
        #Process input z argument to accept ranges or individual numbers for layers
        args.z = args.z.split(',')
        layers = []

        for z in args.z: 
            if '-' in z:
                layers += range(*map(int,z.split('-')))
            elif ':' in z:
                layers += range(*map(int,z.split(':')))
            else: 
                layers += [int(z)]
       

        # Create one output file per desired layers
        for z in layers:

            # limit padding from reaching invalid Z layers
            blue_l  =  max([z-args.blue_padding, 0])
            blue_h  =  min([z+args.blue_padding,32]) 
            green_l =  max([z-args.green_padding,0]) 
            green_h =  min([z+args.green_padding,32]) 

            try:
                r = red[z,:2139,:2139]
                b = blue  [blue_l : blue_h, :2139,:2139].mean(axis = 0) #.astype('uint16') 
                g = green [green_l :green_h,:2139,:2139].sum( axis = 0).astype('uint16') 

            except IndexError: 
                        print(f'ERROR: a requested Z-layer does not exist in one of input images. Constructive image has dimension: {green.shape}')
                        exit()


            # Small blur of constructive signal
            # Limit min and max values (a pseudo-binarization)
            #g[g>10] = 10                      # Before gaussian intended 
            #g = cv2.GaussianBlur(g,(5,5), 0)  
            #g[g<4] = 0                    


            ############# Illumination correction ############
            if args.correct_ilum:     
                thresh = 250 * args.correct_ilum # Arbitrary threshold to avoid areas outside tissue to be normalized up

                # Blurs used for tile normalization
                #rblur = cv2.GaussianBlur(r,(501,501),0)          
                bblur = cv2.GaussianBlur(r,(51,51),0)   

                # Normalize WGA staining based on overall exposition flow (with large gaussian blur)  (DAPI does not work well)
                #rblur [rblur < thresh] = thresh       
                #r = r / rblur

                b = b/bblur    

                # Anyway limit maximum WGA
                #r [r> np.percentile(r,99)] = np.percentile(r,99)

            ##################################################

            # Expand range to 0-255 for max contrast
            r = ((r - r.min()) / (r.max()-r.min()) * 255 ).astype('uint8')
            g = ((g - g.min()) / (g.max()-g.min()) * 255 ).astype('uint8')
            b = ((b - b.min()) / (b.max()-b.min()) * 255 ).astype('uint8')

            # Stack channels and write png 
            rgb = np.dstack([b,g,r])  # stacks 3 h x w arrays -> h x w x 3
            cv2.imwrite(f'{outname}_z{z}_RGB.png',rgb,[int(cv2.IMWRITE_PNG_COMPRESSION), 100])
            
            print(f'Saving: {outname}_z{z}_RGB.png') 


    else:

        g = green 
        r = red 
        b = blue 

        # Small blur of constructive signal
        # Limit min and max values (a pseudo-binarization)
        #g[g>10] = 10                      # Before gaussian intended 
        #g = cv2.GaussianBlur(g,(5,5), 0)  
        #g[g<4] = 0                    

        ##################################################

        # Expand range to 0-255 for max contrast
        r = ((r - r.min()) / (r.max()-r.min()) * 255 ).astype('uint8')
        g = ((g - g.min()) / (g.max()-g.min()) * 255 ).astype('uint8')
        b = ((b - b.min()) / (b.max()-b.min()) * 255 ).astype('uint8')

        # Stack channels and write png 
        rgb = np.dstack([b,g,r])  # stacks 3 h x w arrays -> h x w x 3
        print(f'Saving: {outname}_RGB.png') 
        cv2.imwrite(f'{outname}_RGB.png',rgb,[int(cv2.IMWRITE_PNG_COMPRESSION), 100])
        
