# DSP Assignment 4b
Samuel Stark (sws35)

This project requires the MATLAB Signal Processing Toolbox to function.
I also had the DSP System Toolbox, Image Processing Toolbox, and Symbolic Math Toolbox installed during development,
but I don't think they are needed to run this program.

## Files
- assignment4.m contains the root function for classifying a file, 
    as well as a function for computing compression quality from quantization.
- assignment4_test.m calls assignment4() on the provided test images.
- assignment4_train.m calls assignment4() on the training images I generated.
- /training/ stores the converter.py script and the converted image files I generated.
- figures_fft.m and assignment4_figs.m reproduce figures from the report.
- identify_qs.m is called by assignment4() to find the quantization levels for an image.
- estimate_dc_fft_qs.m is called by identify_qs() to find the quantization levels for an image.
- image_DCT_params.m reads in an image and returns the data for the luminance component.
- neelamani2006_select_best_q.m evaluates the Neelamani2006 grayscale image CHEst algorithm. 
    I wanted to use it in assignment4, but it wasn't very helpful in the end.