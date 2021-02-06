# Open Cell Foam Modeling Using Voronoi Diagram

- Voronoi diagram was generated using n seeds in unit cube to model open-cell foam.
- Voronoi seeds were generated using [Modified MPS](http://www.cs.sandia.gov/~samitch/papers/cccg-present.pdf) for regularity less than 80%.
- Voronoi seeds were generated using [Simple MPS by Implicit Quad-Trees](https://link.springer.com/chapter/10.1007%2F978-3-662-44900-4_13) for regularity higher than 80%.
- For every regularity and relative density, 20 models were generated by a different random distribution of seeds. 
- These combinations generated 300 models.
- The code was developed to generate these models automatically.
- This program encapsulates [TetGen library](http://wias-berlin.de/software/tetgen/) to get Voronoi cells. 
- The code was used to generate virtual structure models which simulate the internal structure of foam material. For the visualization of the simulation [watch this video on Youtube](https://www.youtube.com/watch?v=NMC2FxQ047E).  
![Simulation Output](https://academy.3ds.com/sites/default/files/Graphical%20abstract_0.jpg)

## Compiling
The software has been compiled using cmake and g++ on Linux. The code uses c++17 extensions.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory in the top level directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./foam`. 

## Interface
There is no interface for the current version. You can modify the output by modifing the main function.  

## Copyright and Credit
To credit this work, please cite the following paper and software:

@article{author = {A. M. Fathy; M. H. Abdelshafy; M. R. A Atia},
 title = {MODELING OPEN-CELLED ALUMINUM FOAMS STRUCTURE USING 3-D VORONOI DIAGRAM},
 journal = {AMME},
 volume = {18},
 year = {2018},
 articleno = {76},
 numpages = {12},
 url = [amme](https://amme.journals.ekb.eg/article_35022.html),
 doi = [AMME](10.21608/AMME.2018.35022)
 keywords = {Open-celled Aluminium foam; Computer simulations; Finite element modelling},
   note={Open source software available from \url{https://github.com/ahmedfathy17/Open-cell-Foam-Model-Generation-Using-3-D-Voronoi-diagram.git}}
} 


