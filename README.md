# Julia_imzML_GUI<br />
A Graphical User Interface created with Julia Genie Builder made for this repository: https://github.com/CINVESTAV-LABI/julia_mzML_imzML that generates images with the imzML and ibd files.<br />
For this to work you must insert the next lines as they are depicted in your Julia terminal for the first time only: <br />
```
julia
]
add Pkg Libz ; add https://github.com/CINVESTAV-LABI/julia_mzML_imzML ; add PlotlyBase
```
Then open a terminal on the directory "JuliaIMZML_GUI" and put the next code: <br />
```
julia
]
activate .
##(press backspace)
using Genie
##(an option may appear, select Y to download genie if that happens)
Genie.loadapp() ; up()
##(if it isn't the first time you're running the program:)
using Genie ; Genie.loadapp() ; up()
```
And finally open the port that gets generated to access to the GUI<br />
