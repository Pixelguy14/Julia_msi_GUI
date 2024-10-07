# Julia_msi_GUI<br />
A Graphical User Interface for IMS in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## Local installation
1. Launch Julia
```
julia
```
2. set working directory to Julia_msi_GUI (this repository) code using
```
cd("PathToRepository/Julia_msi_GUI-main")
```
3. Enter Pkg mode by pressing the close square bracket once Julia has been loaded: ** Mostrar imagen del modo pkg**
```
]
```
3. Install the following libraries
```
add Pkg Libz ; add https://github.com/CINVESTAV-LABI/julia_mzML_imzML ; add PlotlyBase ; add Statistics ; add CairoMakie ; add Colors ; add Genie
```

## Load user interface
1. set working directory to Julia_msi_GUI (this repository) code using
```
cd("PathToRepository/Julia_msi_GUI-main")
```
2. In the Pkg mode (enter by pressing "]" once Julia is rining) enter:
```
activate .
```
3. Get out from Pkg mode by pressing the backspace key and enter next line:
```
using Genie ; Genie.loadapp() ; up()
```
6. And finally open the port that gets generated by Genie to access to the GUI (https://localhost:????)<br />
