# High-Resolution-Optical-Mapping
Open and run opticalmapping.m in Matlab. For processing and analysis of high-resolution optical mapping data of the heart. The graphical user interface (GUI) should pop up.
Initially, the option of quickly representing the raw data allows for quick assessment of the signal quality. The user can either choose to process and anayze the entire image or select a region of interest (ROI). 
First step is signal processing which includes baseline drift correction for bleaching. The user can either use the default settings or change the parameters for PCA processing, LPF, MA, and smoothing filters. 
Next, user can generate colourmaps, calculate conduction velocity, create conduction velocity vector field plots, and optical action potentials.
The results including conduction velocity values for each pixel in a user defined region of interest (polygon), APD50, APD70 and APD90 values are outputted in an excel sheet.
