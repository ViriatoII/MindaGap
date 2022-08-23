# MindaGap
   Takes a single panorama image and fills the empty grid lines with neighbour-weighted values.
   Small box sizes yield limited results but work the best with a high loop number (like 40).  Increase boxsize to overcome bigger gaps.
   
   
INSTALLATION 
------------    

 - You need a terminal and Python installed. If you don't, follow the following:
   - I recommend the [git bash](https://github.com/git-for-windows/git/releases/download/v2.37.2.windows.2/Git-2.37.2.2-64-bit.exe) terminal to run python. Install that. 
   - Right click on the folder you want to work on. A new option exists "Git Bash here", click there. A bash terminal is open.
   - Type python and enter. An automatic windows prompt allows you to install it. IMPORTANT: Tick the box that says "Add to PATH"

To install the python dependencies, write the following on the terminal:

```pip3 install tifffile, opencv-python  ```
   
USAGE  
-----------
Open git bash terminal on your desired directory and run:    
 ```python mindagap.py  <PANORAMA.tif> <BOXSIZE> <LOOPNUM> --edges <True|False> ```

    --edges is optional parameter to blur area around grid, for smoother transitions between tiles with different exposures (EXPERIMENTAL)
   
   
Create RGB composite panorama from gapfilled images  
-----------
Use the extra script, like this:

 ```python rgb_from_z_tiles.py  -b <DAPI.tiff> -r <red_channel.tiff> -g <constructive_green_channel.tiff>  ```



    27/06/2022
    Ricardo Guerreiro
    Resolve Biosciences
    
    
