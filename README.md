# JuliaMSI

https://codeberg.org/LabABI/JuliaMSI

## Local Installation
1. Make sure you have Julia installed and at least Julia version 1.6, If not, follow the next installation guide which contains the latest files and instructions for [juliaup](https://github.com/JuliaLang/juliaup) (recommended) alternatively, you can follow the official Julia installation instructions at [julia](https://julialang.org/downloads/platform/).
2. Download and decompress the following repository (You can download the repository from [download](https://codeberg.org/LabABI/JuliaMSI/archive/main.zip) or locate the option manually at the top right of the page, in the box with [•••] (More operations). For a more advanced method using Git, use this link:  https://codeberg.org/LabABI/JuliaMSI.git).
3. Ensure that you have decompressed the ZIP file and know the working directory where it is located (the folder on your computer where you saved the JuliaMSI files).

## Load User Interface
1. Set the working directory to where JuliaMSI is located (the specific location of the folder on your computer where you decompressed the ZIP) in your terminal.<br>
   Use the cd command (which stands for 'change directory') to navigate to the folder where you unzipped the JuliaMSI files.<br>
   Command for Linux:
   ```
   cd (your directory path of this repository)/juliamsi
   ```
   Command for Windows(Powershell or CMD)/Mac:
   ```
   cd (your directory path of this repository)\juliamsi
   ```
   Example of a correct route:<br>
   ~/Downloads/JuliaMSI-main/juliamsi
2. Without entering the Julia environment, launch the project in your terminal with the following command (which works for all operating systems):
   ```
   julia --threads auto --project=. start_MSI_GUI.jl
   ```
3. After the script has finished loading, you can open a [page](http://127.0.0.1:1481/) in your browser with the web app running.

Additional notes:<br>
After the first boot, when the packages are initialized on your computer, subsequent uses of the app should take less time to load.<br>
It is normal to see errors before initialization; you can refresh the page once it opens to see the complete displayed interface.
Minimum system requirements: 4 core processor, 8 GB RAM<br>

JuliaMSI is a Graphical User Interface for a library of MSI tools in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## License

JuliaMSI is published under the terms of the MIT License.

## Example data

Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
Zenodo. <https://doi.org/10.5281/zenodo.10084132>