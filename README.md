# JuliaMSI<br />
A Graphical User Interface for MSI in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## Local Installation
1. Make sure you have Julia installed, if not, follow the next installation guide for [juliaup](https://github.com/JuliaLang/juliaup) (recommended) or [julia](https://julialang.org/downloads/platform/)
2. Download this repository.

## Load User Interface
1. Set working directory to JuliaMSI (this repository) in your terminal using:<br>
   Linux:
   ```
   cd PathToRepository/JuliaMSI
   ```
   Windows/Mac:
   ```
   cd PathToRepository\JuliaMSI
   ```
2. Launch the Julia project in your terminal using:
   ```
   julia --project=. start_MSI_GUI.jl
   ```
3. After the script has finished loading, You can open a page (http://127.0.0.1:1481/) in your browser with the web app running.

Additional notes:<br>
After the first boot initializes the packages in your computer, subsequent uses of the app should take less time to load.<br>
Recomended system requirements: 4 core processor, 8 GB ram<br>

## Example data

Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
Zenodo. <https://doi.org/10.5281/zenodo.10084132>