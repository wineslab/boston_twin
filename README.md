![alt text](images/BOS_G_5_render.png "Boston Twin")

# BostonTwin
Repository for the BostonTwin dataset API.

## Requirements
The basic requirements for the BostonTwin API are based on those of [Sionna]<https://nvlabs.github.io/sionna/>, plus some georeferencing libraries:
1. `python>=3.8`
2. `geopandas`
3. `open3d`
and the corresponding dependencies.\
We provide a requirement file for pip (`requirements.txt`) and conda (`environment.yaml`) to create a python virtual environment with all the dependencies.
Additionally, we provide a DOCKERFILE to generate a container with all the required libraries. This is the preferred method.

## Quickstart
1. Clone this repo
2. Install the requirements. We suggest using 
3. Download the BostonTwin dataset from <http://hdl.handle.net/2047/D20623157> into the `bostontwin` folder.
4. Run the `bostontwin_demo` Jupyter Notebook to see how to use BostonTwin, the Digital Twin of Boston!

## Workflow

![alt text](images/workflow.png "Workflow")

BostonTwin contains the 3D models of the structures [1] and of the antennas [2] in Boston, MA, and relies on [Sionna](<https://nvlabs.github.io/sionna/>) to provide a realistic characterization of the propagation of the electromagnetic signal in the area.

The API offers geo-referencing tools to interact with and manipulate the digital twin. Please refer to the `bostontwin_demo` Jupyter Notebook to see how to use BostonTwin.

[//]: # "## Documentation"
[//]: # "Please refer to the Jupyter notebook"

## Data

![alt text](images/LOD_img.png "LOD")


## Credits
[1] The 3D models of the structures in Boston are published by the Boston Planning and Development Agency (BPDA) [<https://www.bostonplans.org/3d-data-maps/3d-smart-model>].
[2] The antenna location was obtained from the City of Boston's Open Data hub [<https://data.boston.gov/>].
