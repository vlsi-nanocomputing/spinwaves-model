[![DOI](https://zenodo.org/badge/409162306.svg)](https://zenodo.org/badge/latestdoi/409162306)

# Spinwave computational model for all-magnon circuits

The Matlab software offers the possibility to study and design all-magnon circuits taking into account their physical behavior.
This approach enables the architectural exploration of spinwave-based circuits from an architectural point of view.

Different circuits are already available and different parameters can be extracted.

## Getting Started

The entire toolbox is developed on Matlab.

It is necessary to have Matlab installed with the following toolbox:
- Signal Processing Toolbox

### Run the toolbox
The code is structured in different folders:
- GUI
- simulation
- performance_analysis

Within the _GUI_ folder, you can find the file `mainWindow.mlapp`. Run this script to launch the mainwindow of the tool, where you can select the outputs to be printed, select the circuit to study, set the input stimuly and choose between the YIG100 and YIG30 technology node. **Before starting the simulation from the GUI, please move to the root folder of the repository.**

The _simulation_ folder contains the core of the computational model for the two technology nodes. In addition, it is possible to find also the description of the supported circuits.

The _performance_analysis_ folder contains the script to extract the metrics (area, delay, energy).

## Contributing

Branch naming:
- **feature-[FEATURE NAME]**: branch naming convention for new features 
- **fix-[VERSION NUMBER]**: branch naming convention for bug fix related to particular released version
- **fix-[FIX NAME]**: branch naming convention for other bug fix

When a new feature or fix is completed, please submit a pull request.

## Versioning

We use M.m.p as a versioning system. 
- M = Major
- m = minor
- p = patch

There is always back compatibility among minor version and patches only fix possible bugs.
