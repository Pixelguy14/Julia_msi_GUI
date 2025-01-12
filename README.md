# Julia_msi_GUI<br />
A Graphical User Interface for IMS in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## Local Installation
1. Make sure you have Julia installed, if not, follow the next installation guide for [juliaup](https://github.com/JuliaLang/juliaup) (recommended) or [julia] (https://julialang.org/downloads/platform/)
2. Access the repository using this link:<br>
   https://github.com/Pixelguy14/Julia_msi_GUI<br>
   click "<b><> Code</b>" button then hit "<b>Download ZIP</b>"<br>
   unzip the file in your desired location<br>

## Load User Interface
1. Launch Julia in your terminal using:
```
julia
```
2. Set working directory to Julia_msi_GUI (this repository) code using:
```
cd("PathToRepository/Julia_msi_GUI-main")
```
3. With the active directory being Julia_msi_GUI run the next script:
```
include("start_MSI_GUI.jl")
```
4. Alternatively, you can open a terminal in the path of the repository and run the command and skip the next step:
```
julia start_MSI_GUI.jl
```
5. After the script has finished loading, it should open a page in your browser with the web app running.