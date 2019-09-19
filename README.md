[![Build Status](https://travis-ci.org/vrastil/FastSim.svg?branch=master)](https://travis-ci.org/vrastil/FastSim)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/vrastil/FastSim.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/vrastil/FastSim/alerts/)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/vrastil/FastSim.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/vrastil/FastSim/context:cpp)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/vrastil/FastSim.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/vrastil/FastSim/context:python)


# Fast-Sim

## Install
If you are familiar with docker or singularity the easiest way is to use an image with already build program and all dependencies resolved. You can get it from Singularity Hub
````sh
singularity pull shub://vrastil/FastSim-Container
````

### Singularity
TODO: about singularity

TODO: how to build an image

Shell into the container:
````sh
singularity shell shub://vrastil/FastSim-Container
````
Run the container:
````sh
singularity run shub://vrastil/FastSim-Container
````
Build using as a base:
````sh
sudo singularity build FastSim-Container.simg shub://vrastil/FastSim-Container
````
## Usage

### Generate configuration (input)

### Generate initial conditions

### Integration

### Modify gravity solver

### Output files

### Process data

## Documentation
Documentation is available through [github-pages](https://vrastil.github.io/FastSim/).

You can also build it localy with
````sh
cd doc && doxygen Doxyfile
````
You must have a `doxygen` and `graphviz` packages. 

## Results

## License
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.