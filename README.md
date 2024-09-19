# Julia_msi_GUI<br />
A Graphical User Interface created with <i>Julia Genie Builder</i> made for to assist in the use of this repository: https://github.com/CINVESTAV-LABI/julia_mzML_imzML which generates images with the mzMl, imzML and ibd files.<br />
For the correct implementation and use of this repository one must insert the next lines just as they are depicted in your <b>Julia</b> terminal to install the dependencies needed, this part is only needed once. <br />
```
julia
]
add Pkg Libz ; add https://github.com/CINVESTAV-LABI/julia_mzML_imzML ; add PlotlyBase ; add Statistics ; add CairoMakie ; add Colors ; add
add Genie
```
Then open a <b>Julia</b> terminal on the directory <i>"Julia_msi_GUI"</i> and put the next code, this is needed every time you want to run the GUI: <br />
```
]
```
to enter Pkg mode.
```
activate .
```
Press backspace key to exit Pkg mode.
```
using Genie ; Genie.loadapp() ; up()
```
And finally open the port that gets generated to access to the GUI<br />
