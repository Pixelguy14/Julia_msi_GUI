# Julia_imzML_GUI<br />
A graphical User Interface created with Julia Genie Builder made for this repository: https://github.com/CINVESTAV-LABI/julia_mzML_imzML that generates images with the imzML and ibd files.<br />
For this to work you must insert the next lines as they are depicted in your Julia terminal for the first time only.<br />
julia<br />
]<br />
add Pkg Libz<br />
add https://github.com/CINVESTAV-LABI/julia_mzML_imzML<br />
then open a terminal on the directory "JuliaIMZML_GUI"and type: <br />
julia<br />
]<br />
activate .<br />
(press backspace)<br />
using Genie (an option may appear, select Y to download genie if that happens)<br />
Genie.loadapp()<br />
up()<br />
and open the port that gets generated to access to the GUI<br />
