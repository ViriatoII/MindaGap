doc = ''' USAGE:   python duplicate_finder.py  <geneXYZ.txt> <tilesize> 
    (EXPERIMENTAL)

   Takes a single XYZ_coordinates.txt file and searches for duplicates along grid happening at every 2144 pixels.

   Several pixel distance thresholds are still hardcoded.
    - The window of search around the grid is 30 pixels to each side.
    - The maximum perpendicular shift to consider 2 transcripts potential partners is 7 pixels
    - The maximum parallel distance between potential partners is 25 (this is not allowed: O___10____| |____15_____O but this is: O__3_| |____15_____O )
    - The minimum number of transcripts having the mode XYZ shift is 5.


   Steps:
   - It identifies the most common XYZ distance between pairs of potential duplicates across grid (discounting the most common gene)
   - It subtracts the mode XYZ vector from the positions of the second tile
   - If any tile1 transcript is within X distance from a substracted tile2 transcript of same gene, it is flagged as duplicate.
   - The input dataframe is written as output, and duplicates are marked as such in the gene name.
   
   
    02/1/2022
    Ricardo Guerreiro
    Resolve Biosciences'''


def mode3D(xs,ys,zs):
    ''' Returns the most common combination of 3D coordinates'''

    # Transform coords into string for counting
    coords3D = np.array([pair_xdists, pair_ydists, pair_zdists]).T.tolist()
    str_transform =  [str(x)[1:-1] for x in coords3D]
    
    # Get most common coordinate string and return it as numbers
    count = pd.Series(str_transform).value_counts()

    # If there is a matrix with counts
        ## If the most common count is bigger than 5  (maybe bigger than 10 is better?)
    if len(count):
        if count [0] > 5: 

            best = count.index[0].split(',')
            return(np.array([int(c) for c in best]))

    # Reaches this point if almost none or none potential duplicates are found
    print('Probably no duplication between tiles')
    return([0,0,0])

def find_pot_partners (xyzDF, bool_IDs, direction ='horizontal' ):
    ''' Look for partners (same gene) in other side of the grid'''
    partners,dists = {}, {}
    pair_xdists, pair_ydists, pair_zdists = [] ,[], []                                 
    common_gen = xyzDF.gene.value_counts().index[0]

    for i,(og_i,x1,y1,z1,gen1,nada,partners) in xyzDF[bool_IDs].iterrows():

        # Check neighbours before and after (it's sorted by y values)
        for n in [*range(-30,0)]+[*range(1,30)]:  #-1,-2       
            if i +n > 0 and i+n < pNUM and xyzDF.y[i+n]>y1-7:

                if xyzDF.y[i+n]>y1+7: break # Quit searching if neighbors are far in y axis

                partner = xyzDF.loc[i+n]

                # Partner is good if same gene, on the other side of grid, is within 10 pixels in Y direction and 5 in Z.
                if partner.gene == gen1 != common_gen and partner.x > tilesize:  
                    if abs(partner.z - z1) <7: 

                        # Separated handeling for vertical and horizontal grid lines 
                        if direction == 'vertical':
                            if abs(partner.y - y1) <7: 
                                if partner.x - x1  >3 and x1 - partner.x > -25 : # Potential partner has been found  

                                    pair_xdists.append(partner.x - x1) 
                                    pair_ydists.append(partner.y - y1) 
                                    pair_zdists.append(partner.z - z1) 

                                    partners.append(partner['index']) 
                                    xyzDF.loc[old_ids[partner['index']],'partners'].append(og_i) 

                        elif direction == 'horizontal':
                            if abs(partner.x - x1) <7: 
                                if partner.y - y1  >3 and y1 - partner.y > -25 : # Potential partner has been found  

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
    

#################################################################################

                        
if __name__ == '__main__':
    import os,glob, argparse
    import numpy as np, pandas as pd
    import matplotlib.pyplot as plt

   
    parser = argparse.ArgumentParser(description = 'Takes a single XYZ_coordinates.txt file and searches for duplicates along grid happening at every 2144 pixels')
    parser.add_argument("input", help="Input gene xyz coordinates file")
    parser.add_argument("tilesize", nargs = '?', type=int, default = 2144,  help="Tile size (distance between gridlines)")
    args=parser.parse_args()


    # Sanity checks of input
    pathname, extension = os.path.splitext(args.input)

    if os.path.exists(args.input) == False:
        print("Input file does not exist!")
        exit()

    # Read input
    df = pd.read_csv(args.input, sep = '\t', header = None)  # '../../../13_rois_for_demo/panoramas/RiceRoot_Transcripts.txt', sep = '\t', header = None)
    df.columns = ['x','y','z','gene','qual']
    
    print('read input xyz dataframe')

    tilesize = args.tilesize #2144
    panoYmax, panoXmax = df.y.max(), df.x.max()
    duplicated = [] 

    # First iteration going through vertical lines 
    for yjump in range(0, panoYmax, tilesize):
        ymin,ymax = yjump,yjump + tilesize

        for xjump in range(tilesize, panoXmax, tilesize):

            xmin,xmax = xjump - 30, xjump + 30

            # Filtered points for only one grid line 
            df1 = df[(df.x > xmin)*(df.x < xmax) *
                     (df.y > ymin)*(df.y < ymax)].sort_values('y').reset_index()

            pNUM = len(df1)
            left = df1.x < xjump
            df1['partners'] = np.empty((len(df1), 0)).tolist()
   
            old_ids = {k:i for i,k in enumerate(df1['index'])}

            # Find potential partners and report the XYZ distances between them
            pair_xdists,pair_ydists,pair_zdists = find_pot_partners(df1, left, 'vertical')

            # Use original index again and find mode of XYZ shift between potential duplicates
            df1 = df1.set_index('index')
            shift = mode3D(pair_xdists, pair_ydists, pair_zdists)

            # If these is a shift 
            if max(shift) > 2:

                #### Confirm duplications through mutual best partner and distance #####
                for i in df1.index[left]:

                    p1 = pointi(i,shift)
                    bp1 = p1.best_partner
                    
                    if bp1:
                        p2 = pointi(bp1,shift)
                        bp2 = p2.best_partner

                        if i == bp2 and p1.best_multdist < 20:
                            duplicated.append(i)     
            

    # Second iteration going through horizontal lines    
    for xjump in range(0, panoXmax, tilesize):

        xmin,xmax = xjump,xjump + tilesize

        for yjump in range(tilesize, panoYmax, tilesize):

            ymin,ymax = yjump -30, yjump + 30

            # Filtered points for only one grid line 
            df1 = df[(df.x > xmin)*(df.x < xmax) *
                     (df.y > ymin)*(df.y < ymax)].sort_values('y').reset_index()

            pNUM = len(df1)
            bottom = df1.y < yjump
            df1['partners'] = np.empty((len(df1), 0)).tolist()
   
            old_ids = {k:i for i,k in enumerate(df1['index'])}

            # Find potential partners and report the XYZ distances between them
            pair_xdists,pair_ydists,pair_zdists = find_pot_partners(df1, bottom, 'horizontal')

            # Use original index again and find mode of XYZ shift between potential duplicates
            df1 = df1.set_index('index')
            shift = mode3D(pair_xdists, pair_ydists, pair_zdists)
            
            # If these is a shift 
            if max(shift)> 2:
                #### Confirm real duplications by mutual best partner and max distance #####
                for i in df1.index[bottom]:
                    
                    p1 = pointi(i,shift)
                    bp1 = p1.best_partner
                    
                    if bp1:
                        p2 = pointi(bp1,shift)
                        bp2 = p2.best_partner
                        
                        #if p1.point.gene != df1.gene.value_counts().index[0]:

                        if i == bp2 and p1.best_multdist < 20:   # multdist is distance multiplied by weight vector [6,6,1], so Z variation is 6 times less relevant than XY variation
                            duplicated.append(i ) 


                #print(xmin,xmax, ymin,ymax, pNUM, len(duplicated))            
            

    # Replace Gene name by Duplicated and write new XYZ.txt file 
    df.loc[duplicated, 'gene'] = 'Duplicated'
    df.to_csv(pathname + '_markedDups.txt', sep = '\t', header=None, index = False )
