# MindaGap
   Takes a single panorama image and fills the empty grid lines with neighbour-weighted values.
   Small box sizes yield limited results but work the best with a high loop number (like 40).  Increase boxsize to overcome bigger gaps.
 
   Ricardo Guerreiro,   <b> [Resolve Biosciences GmbH](https://resolvebiosciences.com/)   </b>
   
   
   
INSTALLATION 
------------    

 - You need a terminal and Python installed. If you don't, do the following:
   - I recommend the [git bash](https://github.com/git-for-windows/git/releases/download/v2.37.2.windows.2/Git-2.37.2.2-64-bit.exe) terminal. Download and install it.
   - Right click on any folder you want to work on. A new option exists "Git Bash here", click there. A bash terminal is open.
   - Type python and enter. An automatic windows prompt allows you to install it.
   - Enter this in the bash terminal: ```echo "alias python=' winpty python.exe'" >> ~/.bashrc ; source ~/.bashrc```

- To install the python dependencies, enter the following on the terminal:

```pip3 install tifffile opencv-python  ```

- Install MindaGap:

``` cd ; mkdir Programs; cd Programs ```     
``` git clone https://github.com/ViriatoII/MindaGap/  ```


   
USAGE  
-----------
Open git bash terminal on your desired directory and run:    
 ```python ~/Programs/MindaGap/mindagap.py  INPUT_PANORAMA.tif  [BOXSIZE] [LOOPNUM] [--edges True/False] ```

      Optional parameters:
      BOXSIZE Default 3. A larger number allows to overcome large gaps, but makes looses fine details in new filled grid.      
      LOOPNUM Default 40. A smaller number is faster, but the result is less good.       
      --edges Default False. An optional parameter to blur area around grid, for smoother transitions between tiles with different exposures (EXPERIMENTAL).   
   
   
 # Additional scripts
   
### Create RGB composite panorama from gapfilled images  

 ```python ~/Programs/MindaGap/rgb_from_z_tiles.py  -b DAPI.tiff -r red_channel.tiff -g constructive_green_channel.tiff  ```

### Mark duplicate reads along gridline edges 

 ```python ~/Programs/MindaGap/duplicate_finder.py  XYZ_coordinates.csv [Xtilesize] [Ytilesize] [windowsize] [maxfreq] [minMode] [-p True/False]   ```



    
    

Docker  
-----------
To use MindaGap with Docker, you can pull the prebuild docker image from Dockerhub or follow the steps below to build the Docker image.

You can then execute the script using a docker container with all required dependencies available:

```
## Pull the prebuild docker image from Docker hub
docker pull rguerr/mindagap:latest
```

```
#### To build the docker image yourself
## Clone this git repository
git clone https://github.com/ViriatoII/MindaGap.git
cd MindaGap

## Build docker image
docker build -t mindagap .

## Test that docker build worked successfully
docker run --rm -it mindagap:latest \
python /mindagap/mindagap.py
```


How to cite 
--------------
      Ricardo Guerreiro, Florian Wuennemann, & pvtodorov. (2023). ViriatoII/MindaGap: v0.0.3 (0.0.3). Zenodo. https://doi.org/10.5281/zenodo.8120559
