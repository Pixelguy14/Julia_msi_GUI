# JuliaMSI

https://codeberg.org/LabABI/JuliaMSI

## Local Installation
1. Make sure you have Julia installed. JuliaMSI supports **Julia 1.11 and 1.12**. If Julia is not installed, follow the [juliaup](https://github.com/JuliaLang/juliaup) guide (recommended), or the official Julia installation instructions at [julia](https://julialang.org/downloads/platform/).
2. Download and decompress the following repository (You can download the repository from [download](https://codeberg.org/LabABI/JuliaMSI/archive/main.zip) or locate the option manually at the top right of the page, in the box with [•••] (More operations). For a more advanced method using Git, use this link:  https://codeberg.org/LabABI/JuliaMSI.git).
3. Ensure that you have decompressed the ZIP file and know the working directory where it is located (the folder on your computer where you saved the JuliaMSI files).
4. Julia version 1.11 or 1.12 is required.
5. Ensure that the project's `public/` directory is writable. JuliaMSI writes generated images, masks, and `registry.json` there and will stop with a clear error if it cannot.

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
2. Install/verify the Julia dependencies once (safe to repeat):
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
3. Without entering the Julia environment, launch the project in your terminal with the following command (which works for all operating systems):
   ```bash
   julia --threads auto --project=. start_MSI_GUI.jl
   ```
4. After the script has finished loading, you can open a [page](http://127.0.0.1:1481/) in your browser with the web app running.

   The GUI starts in production mode by default. For live-reload development use:
   ```bash
   JULIA_MSI_DEV=1 julia --threads auto --project=. start_MSI_GUI.jl
   ```
   Set `JULIA_MSI_VERBOSE=1` to see the otherwise-filtered HTML parser warnings.

Additional notes:<br>
After the first boot, when the packages are initialized on your computer, subsequent uses of the app should take less time to load.<br>
It is normal to see errors before initialization; you can refresh the page once it opens to see the complete displayed interface.
Minimum system requirements: 4 core processor, 8 GB RAM<br>

JuliaMSI is a Graphical User Interface for a library of MSI tools in Julia: https://github.com/CINVESTAV-LABI/julia_mzML_imzML

## Build the System Image (One-time)
Alternatively, you can generate custom pre-compiled `.so` / `.dll` system images by running the build script in your directory. This bakes the heavy mass-spectrometry kernels and UI dependencies into a single machine-code binary for ultra-fast startup times.

### 1. Build the GUI Sysimage (for full App usage)
```bash
julia --project=. build_sysimage.jl
```
This may take 5–15 minutes. The resulting `sys_msi_gui.so` (or `.dll` on Windows) file will be around **300MB–600MB** because it contains the pre-compiled machine code for your entire graphical environment.
*Note: A sysimage built on Linux (.so) will not work on Windows. You must run the build script once on each target operating system.*

Once `sys_msi_gui.so` is created in your directory, adapt your launch command:
```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
# Adjust the .so extension to .dll if on Windows or .dylib if on macOS
julia --threads auto --project=. --sysimage sys_msi_gui.so start_MSI_GUI.jl
```

### 2. Build the Headless Sysimage (for Scripts & Data Scientists)
If you only want to run scripts or the core `MSI_src` engine without the overhead of the Genie web server and Plotly, you can build a stripped-down, high-speed system image:

```bash
julia --project=. build_sysimage.jl --headless
```

Launch your automated tests or custom processing scripts using the headless image to bypass virtually all compilation overhead:
```bash
julia --threads auto --project=. --sysimage sys_msi_headless.so test/test_streaming_pipeline.jl
```

### 3. Docker Option for High-Speed Processing
For maximum portability and execution speed without clutter, we provide a `Dockerfile.headless`. This option creates a minimal, ultra-fast container that completely omits `app.jl` and the web interface. 

It automatically builds and bundles the headless sysimage, giving researchers an instant-start engine ready for heavy-duty `.imzML` batch pipelines with zero JIT latency.

To build and run the Docker image:
```bash
docker build -t juliamsi-headless -f Dockerfile.headless .
# Run a script mapped from your local data folder:
docker run -it -v /my/local/data:/data juliamsi-headless julia -J sys_msi_headless.so /data/my_processing_script.jl
```

## License

JuliaMSI is published under the terms of the MIT License.

## Example data

Robert Winkler. (2023). mzML mass spectrometry and imzML mass spectrometry imaging test data [Data set]. 
Zenodo. <https://doi.org/10.5281/zenodo.10084132>

## Plotting API note

`plot_slice`, `plot_slices`, and `generate_colorbar_image` are GUI/plotting helpers.
They now live in the GUI layer (`julia_imzML_visual.jl`) so that `using MSI_src`
does not load CairoMakie for headless/console use. They are not exported from
`MSI_src`; calling `MSI_src.plot_slice(...)` returns an actionable migration
error instead of silently failing with `UndefVarError`.

## User Guides

For a general guide on how to use the JuliaMSI platform, please refer to the [JuliaMSI Main Guide](JuliaMSI_MainGuide.md).

For a more detailed guide on the preprocessing pipeline, please refer to the [JuliaMSI Preprocessing Guide](JuliaMSI_Preprocessing.md).
