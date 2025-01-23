# Julia_msi_GUI<br />
A Graphical User Interface for MSI in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## Local Installation
1. Make sure you have Julia installed, if not, follow the next installation guide for [juliaup](https://github.com/JuliaLang/juliaup) (recommended) or [julia](https://julialang.org/downloads/platform/)
2. Access the repository using this link:<br>
   https://github.com/Pixelguy14/Julia_msi_GUI<br>
   click "<b><> Code</b>" button then hit "<b>Download ZIP</b>"<br>
   unzip the file in your desired location<br>

## Load User Interface
1. Set working directory to Julia_msi_GUI (this repository) in your terminal using:
   linux:
   ```
   cd PathToRepository/Julia_msi_GUI-main
   ```
   windows/mac:
   ```
   cd PathToRepository\Julia_msi_GUI-main
   ```
2. Launch the Julia project in your terminal using:
   ```
   julia --project=. start_MSI_GUI.jl
   ```
3. After the script has finished loading, it should open a page (http://127.0.0.1:1481/) in your browser with the web app running.

Additional notes:
After the first boot initializes the packages in your computer, subsequent uses of the app should not take longer to load.
Recomended system requirements: 4 core processor, 8 GB ram
Minimum system requirements: 2 core processor, 8 GB ram (long loading times)