![alt text](images/BOS_G_5_render.png "Boston Twin")

# BostonTwin
Repository for the BostonTwin dataset API.

## Quickstart
1. Clone this repo
2. Download the BostonTwin dataset from <http://hdl.handle.net/2047/D20623157> into the `bostontwin` folder.
3. Run the `bostontwin_demo` Jupyter Notebook to see how to use BostonTwin, the Digital Twin of Boston!

## Workflow

![alt text](images/workflow.png "Workflow")

BostonTwin contains the 3D models of the structures [1] and of the antennas [2] in Boston, MA, and relies on [Sionna]<https://nvlabs.github.io/sionna/> to provide a realistic characterization of the propagation of the electromagnetic signal in the area.

The API offers a set of geo-referencing tools to interact with and manipulate the digital twin. Please refer to the `bostontwin_demo` Jupyter Notebook to see how to use BostonTwin.

## Documentation

## Data

![alt text](images/LOD_img.png "LOD")


## Credits
[1] The 3D models of the structures in Boston are published by the Boston Planning and Development Agency (BPDA) [<https://www.bostonplans.org/3d-data-maps/3d-smart-model>].
[2] The antenna location was obtained from the City of Boston's Open Data hub [<https://data.boston.gov/>].
