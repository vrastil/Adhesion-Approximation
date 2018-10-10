[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1667)

#  FastSim Container
This repository serves as a portble solution to install [FastSim](https://github.com/vrastil/Adhesion-Approximation)

Use singularity recipes in subdirectories to build your own containers

````
sudo singularity build fastsim.simg Singularity
````

or use already build ones at SingularityHub

````
singularity pull --name fastsim.img shub://vrastil/FastSim-Container
````
