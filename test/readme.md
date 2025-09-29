to run the tests:

time julia --project=. test/run_tests.jl

to enable the different test cases, there are boolean variables

test1 = true
test2 = true 
test3 = true

set them to true to perform that specific test.

test 1 consists in mzml validation and spectrum plotting

test 2 consists in imzml conversion and validation

test 3 consists in imzml validation and spectrum and slice plotting

use the global variables:
TEST_MZML_FILE, SPECTRUM_TO_PLOT, CONVERSION_SOURCE_MZML, CONVERSION_SYNC_FILE, CONVERSION_TARGET_IMZML, TEST_IMZML_FILE, MZ_VALUE_FOR_SLICE, MZ_TOLERANCE, COORDS_TO_PLOT, RESULTS_DIR 
to modify input that the test script will recieve.
