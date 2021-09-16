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

Within the GUI folder, you can find the file `mainWindow.mlapp`. Run this script to launch the mainwindow of the tool, where you can select the outputs to be printed, select the circuit to study, set the input stimuly and choose between the YIG100 and YIG30 technology node.

## 1) simulation
This folder contains the developed simulation tool of this thesis.


## 2) performance_analysis
This folder contains all the codes of performance analysis (area, delay, energy).


## 3) synopsys
This folder contains all the files for synopsys synthesis.


## 4) Other_Scripts
This folder contains some useful scripts. 


## 5) guide.pdf
This pdf can be considered as a guide of simulation tool and performance analysis.


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
