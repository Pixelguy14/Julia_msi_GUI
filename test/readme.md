# JuliaMSI Test Suite

This document provides instructions on how to set up and run the test suite for the `JuliaMSI` package. The tests validate the core functionality of the data processing workflows, including loading, converting, and visualizing mass spectrometry data.

## Local Installation and Setup

### 1. Prerequisites
Make sure you have Julia installed (at least version 1.6). If not, we recommend using `juliaup` for installation.
- **juliaup**: [https://github.com/JuliaLang/juliaup](https://github.com/JuliaLang/juliaup)
- **Official binaries**: [https://julialang.org/downloads/](https://julialang.org/downloads/)

### 2. Download the Project
Download and decompress the repository. You can get it from the official Codeberg page or clone it using Git.
- **Download ZIP**: [https://codeberg.org/LabABI/JuliaMSI/archive/main.zip](https://codeberg.org/LabABI/JuliaMSI/archive/main.zip)
- **Git URL**: `https://codeberg.org/LabABI/JuliaMSI.git`

Ensure you know the directory where the `JuliaMSI` files are located.

### 3. Download Test Data
The test script requires example mass spectrometry data files. You can download the required test data from Zenodo:
- **Example Data**: [https://doi.org/10.5281/zenodo.10084132](https://doi.org/10.5281/zenodo.10084132)

Download the data and place it in a known location on your computer.

### 4. Configure the Test Script
Before running the tests, you must edit the `test/run_tests.jl` file to point to the test data you downloaded.

1.  Open `test/run_tests.jl` in a text editor.
2.  Locate the `CONFIG: PLEASE FILL IN YOUR FILE PATHS HERE` section.
3.  Update the constant variables (e.g., `TEST_MZML_FILE`, `CONVERSION_SOURCE_MZML`, `CONVERSION_SYNC_FILE`) with the absolute paths to the corresponding files on your system.

## Running the Tests

1.  **Navigate to the Project Directory**:
    Open your terminal and use the `cd` command to navigate to the root folder of the `JuliaMSI` project.
    ```bash
    cd /path/to/your/JuliaMSI
    ```

2.  **Execute the Test Script**:
    Run the following command from the project's root directory. This will install the necessary dependencies and run the tests.
    ```bash
    julia --project=. test/run_tests.jl
    ```

3.  **Check the Results**:
    The script will print its progress to the console. Any generated images (plots and image slices) will be saved in the `test/results/` directory.

## Test Case Configuration

You can customize the test run by editing the variables in `test/run_tests.jl`.

### Enabling and Disabling Test Cases
You can run or skip specific test cases by setting the corresponding boolean variables to `true` or `false`.

```julia
test1 = true  # Runs Test Case 1
test2 = true  # Runs Test Case 2
test3 = true  # Runs Test Case 3
```

-   **Test Case 1**: Validates a standard `.mzML` file and generates plots for a single spectrum, the total spectrum, and the average spectrum.
-   **Test Case 2**: Tests the conversion of a `.mzML` file (and its corresponding `.txt` sync file) into the `.imzML` format. It then validates the newly created file.
-   **Test Case 3**: Validates an existing `.imzML` file. It generates plots for a spectrum at specific coordinates, the total spectrum, the average spectrum, and an ion image slice.

### Input Configuration Variables
All configuration variables are located in the `CONFIG` section of `test/run_tests.jl`.

-   `TEST_MZML_FILE`: Path to the standard `.mzML` file for Test Case 1.
-   `SPECTRUM_TO_PLOT`: The index of the spectrum to plot from the `.mzML` file in Test Case 1.
-   `CONVERSION_SOURCE_MZML`: Path to the imaging `.mzML` file to be converted in Test Case 2.
-   `CONVERSION_SYNC_FILE`: Path to the `.txt` synchronization file corresponding to `CONVERSION_SOURCE_MZML`.
-   `CONVERSION_TARGET_IMZML`: The output path for the `.imzML` file generated in Test Case 2. This is also used as the default input for Test Case 3.
-   `TEST_IMZML_FILE`: Path to the `.imzML` file to be tested in Test Case 3.
-   `MZ_VALUE_FOR_SLICE`: The m/z value to use for creating an ion image slice in Test Case 3.
-   `MZ_TOLERANCE`: The tolerance (+/-) for the `MZ_VALUE_FOR_SLICE` when generating the image.
-   `COORDS_TO_PLOT`: An `(X, Y)` tuple specifying the coordinates of the spectrum to plot from the `.imzML` file in Test Case 3.
-   `RESULTS_DIR`: The directory where output images will be saved.