doc = ''' USAGE:   python duplicate_finder.py  <geneXYZ.txt> <Xtilesize> [Ytilesize==Xtilesize] [windowsize] [maxfreq] [minMode]  [-p True/False] 
    (EXPERIMENTAL)

   Takes a single XYZ_coordinates.txt file and searches for duplicates along grid happening at every +/- 2144 pixels.

positional arguments:
  input                 Input gene xyz coordinates file
  Xtilesize             Tile size (distance between gridlines) on X axis. Default = 2144
  Ytilesize             Tile size (distance between gridlines) on Y axis. Default = Xtilesize
  windowsize            Window arround gridlines to search for duplicates. Default = 30 (to each side of gridline)
  maxfreq               Maximum transcript count to calculate X/Y shifts (better to discard very common genes here). Default = 400
  minMode               Minumum occurances of ~XYZ_shift to consider it valid. Default = 10 

optional arguments:
  -p PLOT, --plot PLOT  Illustrative lineplot of duplicated pairs with annotated XYZ shift per tileOvlap

   Several pixel distance thresholds are tweekable through user input, but are default as:
    - The windowsize of search around the gridbreak is 30 pixels to each side.
    - The maximum parallel distance between potential partners is 25 (windowsize - 5) (this is not allowed: O___10____| |____15_____O but this is: O__3_| |____15_____O )
    - The maximum perpendicular shift to consider 2 transcripts potential partners is 7 pixels (hardcoded)
    - The minimum number of transcripts having ~ mode XYZ shift is 10.


   Steps:
   - It identifies the most common XYZ distance between pairs of potential duplicates across grid (discounting the most common gene)
   - It subtracts the mode XYZ vector from the positions of the second tile
   - If any tile1 transcript is within X distance from a substracted tile2 transcript of same gene, it is flagged as duplicate.
   - The input dataframe is written as output, and duplicates are marked as such in the gene name.
   
   
    02/11/2022
    Ricardo Guerreiro
    Resolve Biosciences'''


def mode3D(xs,ys,zs, minMode):
    ''' Returns the most common combination of 3D coordinates
    Takes into consideration if 2nd mode and 3rd mode are very close to 1st mode'''

    # Transform coords into string for counting
    coords3D = np.array([pair_xdists, pair_ydists, pair_zdists]).T.tolist()
    str_transform =  [str(x)[1:-1] for x in coords3D]
    
    # Get most common coordinate string and return it as numbers
    count = pd.Series(str_transform, dtype=str).value_counts()

    # If there is a matrix with counts
    if len(count):

        check_rows = min([3, len(count)])
        best = [np.array(x.split(','), dtype=int) for x in count.index[:check_rows]]

        # Sum counts of second and third mode if they are very close to first mode
        for i in range(1,check_rows):
            if abs(best[0] - best[i] ).sum()<3:
                count [0] += count[i] 

        ## If the most common count is bigger than minMode, shift is accepted as valid
        if count [0] > minMode: 

            #print(f'At least {count [0]} duplications between tiles Y{xmin}:{xmax}, X{ymin}:{ymax}')

            best = count.index[0].split(',')
            return(np.array([int(c) for c in best]))

    # Reaches this point if almost none or no potential duplicates are found
    #print(f'Probably no duplication between tiles Y{xmin}:{xmax}, X{ymin}:{ymax}')
    return([0,0,0])

def find_pot_partners_horizontal (xyzDF, bool_IDs,max_freq):
    ''' Look for partners (same gene) in other side of the grid'''
    partners,dists = {}, {}
    pair_xdists, pair_ydists, pair_zdists = [] ,[], []    
    pNUM = len(xyzDF) 

    # Deal with empty input and filter out too common genes
    if len(bool_IDs) > 0: 
        gen_counts = xyzDF.gene.value_counts()

        common_gens = []

        for g, c in gen_counts.items():
            if c > max_freq:
                common_gens.append(g) 
            else: break
    else:  return([0,0,0])
                           

    #print('horizontal') ; #print(xyzDF[bool_IDs].gene.value_counts())

    for i,(og_i,x1,y1,z1,gen1,partners) in xyzDF[bool_IDs].iterrows():

        # Check neighbours before and after (it's sorted by y values)
        for n in [*range(-30,0)]+[*range(1,30)]:  #-1,-2       
            if i +n > 0 and i+n < pNUM and xyzDF.x[i+n]>x1-7:
                if xyzDF.x[i+n]>x1+7: break # Quit searching if neighbors are far in y axis

                partner = xyzDF.loc[i+n]

                # Partner is good if same gene, on the other side of grid, is within 10 pixels in Y direction and 5 in Z.
                if partner.gene == gen1 not in common_gens: # and partner.y >= yjump: 
                    if abs(partner.z - z1) <7 and abs(partner.x - x1) <7: 
                        if partner.y - y1  >2 and partner.y - y1 < w-5: # Do not consider very large parallel distances (far left transcript only connects with transcripts on the left of panel2)

                            # Potential partner has been found  
                            pair_xdists.append(partner.x - x1) 
                            pair_ydists.append(partner.y - y1) 
                            pair_zdists.append(partner.z - z1) 

                            partners.append(partner['index']) 
                            xyzDF.loc[old_ids[partner['index']],'partners'].append(og_i) 

    return(pair_xdists,pair_ydists,pair_zdists)

def find_pot_partners_vertical (xyzDF, bool_IDs,max_freq):
    ''' Look for partners (same gene) in other side of the grid'''
    partners,dists = {}, {}
    pair_xdists, pair_ydists, pair_zdists = [] ,[], []  
    pNUM = len(xyzDF) 
    

    # Deal with empty input and filter out too common genes
    if len(bool_IDs) > 0: 
        gen_counts = xyzDF.gene.value_counts()

        common_gens = []

        for g, c in gen_counts.items():
            if c > max_freq:
                common_gens.append(g) 
            else: break
    else: return([0,0,0])
                        

    for i,(og_i,x1,y1,z1,gen1,partners) in xyzDF[bool_IDs].iterrows():

        # Check neighbours before and after (it's sorted by y values)
        for n in [*range(-30,0)]+[*range(1,30)]:  #-1,-2       
            if i +n > 0 and i+n < pNUM and xyzDF.y[i+n]>y1-7:

                if xyzDF.y[i+n]>y1+7: break # Quit searching if neighbors are far in y axis

                partner = xyzDF.loc[i+n]

                # Partner is good if same gene, on the other side of grid, is within 10 pixels in Y direction and 5 in Z.
                if partner.gene == gen1 not in common_gens and partner.x >= xjump: 
                    if abs(partner.z - z1) <7 and abs(partner.y - y1) <7: 

                        if partner.x - x1  >2 and partner.x - x1 < w -5 : # Do not consider very large parallel distances (far left transcript only connects with transcripts on the left of panel2)

                            # Potential partner has been found  
                            pair_xdists.append(partner.x - x1) 
                            pair_ydists.append(partner.y - y1) 
                            pair_zdists.append(partner.z - z1) 

                            partners.append(partner['index']) 
                            xyzDF.loc[old_ids[partner['index']],'partners'].append(og_i) 

    return(pair_xdists,pair_ydists,pair_zdists)

class pointi:
    ''' Point class with computed distances to partners ''' 
    def __init__(self,idx,shift):
        self.id = idx
        self.point = df1.loc[idx]
        self.partners = self.point.partners
        self.best_partner = None
        
        if len(self.partners)== 0: return(None)
    
        # Signal to deal with p1 and p2 being either first-second or second-first 
        sig = (self.point.values[:3] - df1.loc[self.partners[0]].values[:3]).sum() +0.0001
        sig /= abs(sig)
        
        # Calculate distance between point1 and partners
        self.dists = [abs(sig * self.point.values[:3] - sig*df1.loc[p].values[:3] - shift).astype(int) for p in self.partners]   
        self.dists3D = [d.sum() for d in self.dists]  
        self.best_dist =  min(self.dists3D)
        
        # Identify which partner has best (minimum) distance         
            # Giving more weight to x,y coordinates than z 
        mult_dists = (self.dists * np.array([6,6,1])).sum(axis=1)
        self.best_multdist =    mult_dists.min()
        minNUM = (mult_dists == mult_dists.min()).sum()  # (self.dists3D == self.best_dist).sum() 

        if minNUM == 1:
            self.best_partner = self.partners[mult_dists.argmin()] 
        else:     # Deal with multiple best partners (equal dist)
            
            # Prefer distance arrays with smaller SD ([3,1] is better than [4,0]). 
            min_sd = np.std(self.dists,axis=1).argmin()    # If SD is exactly equal, first is taken
            self.best_partner = self.partners[min_sd]
    

###########################   MAIN  ##############################################

                        
if __name__ == '__main__':
    import os,glob, argparse
    import numpy as np, pandas as pd
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description = 'Takes a single XYZ_coordinates.txt file and searches for duplicates along grid happening at every 2144 pixels')
    parser.add_argument("input", help="Input gene xyz coordinates file")
    parser.add_argument("Xtilesize", nargs = '?', type=int, default = 2144,  help="Tile size (distance between gridlines) on X axis")
    parser.add_argument("Ytilesize", nargs = '?', type=int, default = None,  help="Tile size (distance between gridlines) on Y axis")
    parser.add_argument("windowsize", nargs = '?', type=int, default = 30,  help="Window arround gridlines to search for duplicates")
    parser.add_argument("maxfreq", nargs = '?', type=int, default = 400,  help="Maximum transcript count to calculate X/Y shifts (better to discard very common genes)")
    parser.add_argument("minMode", nargs = '?', type=int, default = 10,  help="Minumum occurances of ~XYZ_shift to consider it valid")
    parser.add_argument("-p", "--plot", default = None,  help="Illustrative lineplot of duplicated pairs with annotated XYZ shift per tileOvlap")
    args=parser.parse_args()

    # Sanity checks of input
    pathname, extension = os.path.splitext(args.input)

    if os.path.exists(args.input) == False:
        print("Input file does not exist!"); exit()

    # Read input
    df = pd.read_csv(args.input, sep = '\t', header = None)  # '../../../13_rois_for_demo/panoramas/RiceRoot_Transcripts.txt', sep = '\t', header = None)
    df = df.drop(df.columns[4:], axis = 1)
    df.columns = ['x','y','z','gene'] #,'qual']
    print('read input xyz dataframe')

    # Prepare parameters 
    w = args.windowsize
    minMode = args.minMode
    max_freq = args.maxfreq
    Xtilesize = args.Xtilesize #2144
    Ytilesize = args.Xtilesize if args.Ytilesize == None else args.Ytilesize 
    panoYmax, panoXmax = df.y.max(), df.x.max()
    duplicated = [] 
    tileOvlaps, totalDups, tilePairs = 0, 0, 0

    # First iteration going through vertical lines 
    for yjump in range(0, panoYmax, Ytilesize):
        ymin,ymax = yjump,yjump + Ytilesize

        for xjump in range(Xtilesize, panoXmax, Xtilesize):

            xmin,xmax = xjump - w, xjump + w

            # Filtered points for only one grid line 
            df1 = df[(df.x > xmin)*(df.x < xmax) *
                     (df.y > ymin)*(df.y < ymax)].sort_values('y').reset_index()

            # Skip rest if no transcripts exist 
            if df1.empty:
                print(f'Found no transcripts within Y{xmin}:{xmax}, X{ymin}:{ymax}')
                continue     

            if args.plot: # Plot gridlines if asked
                rectangle = plt.Rectangle((xmin,ymin), 2*w , ymax -ymin, fill=False, ec="grey", linewidth = 0.3)
                plt.gca().add_patch(rectangle)

            tilePairs +=1 
            left = df1.x < xjump
            df1['partners'] = np.empty((len(df1), 0)).tolist()
   
            old_ids = {k:i for i,k in enumerate(df1['index'])}

            # Find potential partners and report the XYZ distances between them
            pair_xdists,pair_ydists,pair_zdists = find_pot_partners_vertical(df1, left, max_freq)

            # Use original index again and find mode of XYZ shift between potential duplicates
            df1 = df1.set_index('index')
            shift = mode3D(pair_xdists, pair_ydists, pair_zdists, minMode)

            if args.plot: # Annotate plot with XYZ shift between tiles
                plt.text(xmin-400, ymin+1000, str(shift), size = 4)

            # If these is a shift 
            if max(shift) > 2:
                dups = 0 

                #### Confirm duplications through mutual best partner and distance #####
                for i in df1.index[left]:

                    p1 = pointi(i,shift)
                    bp1 = p1.best_partner
                    
                    if bp1:
                        p2 = pointi(bp1,shift)
                        bp2 = p2.best_partner

                        if i == bp2 and p1.best_multdist < 20:
                            duplicated.append(i) 
                            dups +=1  

                            if args.plot:
                                plt.plot([p1.point.x, p2.point.x],
                                    [p1.point.y, p2.point.y], c = 'red', linewidth = 0.2 )

                # Count Tileoverlap with at least 5 duplicated points    
                if dups > 4: 
                    totalDups += dups
                    tileOvlaps +=1 
            

    # Second iteration going through horizontal lines    
    for xjump in range(0, panoXmax, Xtilesize):

        xmin,xmax = xjump,xjump + Xtilesize

        for yjump in range(Ytilesize, panoYmax, Ytilesize):

            ymin,ymax = yjump -w, yjump + w

            # Filtered points for only one grid line 
            df1 = df[(df.x > xmin)*(df.x < xmax) *
                     (df.y > ymin)*(df.y < ymax)].sort_values('x').reset_index()

            if df1.empty:
                print(f'Found no transcripts within Y{xmin}:{xmax}, X{ymin}:{ymax}')
                continue  

            if args.plot: # Plot gridlines if asked
                rectangle = plt.Rectangle((xmin,ymin), xmax -xmin, 2*w, fill=False, ec="grey", linewidth = 0.3)
                plt.gca().add_patch(rectangle)

            tilePairs +=1 
            bottom = df1.y < yjump
            df1['partners'] = np.empty((len(df1), 0)).tolist()
   
            old_ids = {k:i for i,k in enumerate(df1['index'])}

            # Find potential partners and report the XYZ distances between them
            pair_xdists,pair_ydists,pair_zdists = find_pot_partners_horizontal(df1, bottom,max_freq)

            # Use original index again and find mode of XYZ shift between potential duplicates
            df1 = df1.set_index('index')
            shift = mode3D(pair_xdists, pair_ydists, pair_zdists, minMode)

            if args.plot: # Annotate plot with XYZ shift between tiles
                plt.text(xmin+1000, ymin-100, str(shift), size = 4)
            
            # If there is a shift 
            if max(shift)> 2:
                dups = 0 

                #### Confirm real duplications by mutual best partner and max distance #####
                for i in df1.index[bottom]:
                    
                    p1 = pointi(i,shift)
                    bp1 = p1.best_partner
                    
                    if bp1:
                        p2 = pointi(bp1,shift)
                        bp2 = p2.best_partner
                        
                        if i == bp2 and p1.best_multdist < 20:   # multdist is distance multiplied by weight vector [6,6,1], so Z variation is 6 times less relevant than XY variation
                            dups  +=1 
                            duplicated.append(i)

                            # Plot line between point pair if asked for
                            if args.plot: 
                                plt.plot([p1.point.x, p2.point.x],
                                    [p1.point.y, p2.point.y], c = 'red', linewidth = 0.2 )

                # Count Tiles as overlapping if least 10 duplicated points 
                if dups > 9: 
                    totalDups += dups
                    tileOvlaps +=1 
                    #print(xmin,xmax, ymin,ymax, pNUM, len(duplicated))            

    if args.plot: 
        plt.xlabel('X') #;  plt.xlim([df.x.min(), df.x.max()])
        plt.ylabel('Y') #;  plt.ylim([df.y.min(), df.y.max()])


        # Plot path of microscope if -p is an existing file (positions.log file)
        if os.path.exists(args.plot):
            sample_name = args.input.split('_')[-2]

            pos_log = pd.read_csv(args.plot, header = None)
            pos_log.columns = ['time', 'name', 'x','y','z', 'x2','y2','z2']
            pos_log = pos_log[[sample_name in r for r in pos_log.name]]
            pos_log = pos_log[['X_R7' in r for r in pos_log.name]]

            pos_log = pos_log.iloc[:,2:4].values
            dif = pos_log[1][0] -pos_log[0][0]

            l = (pos_log - pos_log.min(axis = 0)) 
            l = l * Xtilesize / dif  + Xtilesize/2


            for i,(x,y) in enumerate(l):
                plt.annotate(i, (x-200,y-100), color = '#3776ab')
                    
                if i< len(l)-1:
                    dX = l[i+1][0] - x 
                    dY = l[i+1][1] - y 

                    plt.arrow(x + dX*0.05,y +dY*0.05 , dX*0.75, dY*0.75, 
                              head_width=100,# width = 0.8,
                              ls='-.', color = '#89cff0') 


        plt.gca().invert_yaxis()
        plt.savefig('XYZshift_visualization.png', dpi = 550)

    # Replace Gene name by Duplicated and write new XYZ.txt file 
    df.loc[duplicated, 'gene'] = 'Duplicated'
    df.to_csv(pathname + '_markedDups.txt', sep = '\t', header=None, index = False )


    # Counts for analysis: Number of tile pairs, number of tile overlaps, total duplicated transcripts
    print('TilePairs, Tile Overlaps, Total duplicated transcripts')
    print(tilePairs, tileOvlaps, totalDups)


