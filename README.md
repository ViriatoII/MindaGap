# MindaGap
   Takes a single panorama image and fills the empty grid lines with neighbour-weighted values.
   Small box sizes yield limited results but work the best with a high loop number (like 40).  Increase boxsize to overcome bigger gaps. 
   
USAGE:   ```python mindagap.py  <PANORAMA.tif> <BOXSIZE> <LOOPNUM> --edges <True|False> ```

    --edges is optional parameter to blur area around grid, for smoother transitions between tiles with different exposures (EXPERIMENTAL)
   
    27/06/2022
    Ricardo Guerreiro
    Resolve Biosciences
