### Files

#### Main Scripts

- `girf_v2.m` contains code to generate a GIRF measurement sequence for use with a Skope field camera.  This mirrors the sequence that was used in our work.

- `generate_GIRF.ipynb` contains a python notebook that demonstrates loading the GIRF measurement data and solving for the GIRF.


#### Other code

- `girf_utils.py` contains helpe functions and the main function used to generate a GIRF from the measurements acquired from the field camera measurements.

- `girf3d_2.m` is a Matlab script for generating a Pulseq sequence to measure a 3D GIRF using the thin slice method.  This was not thoroughly tested (because we switched to the field camera), however initial tests showed that it performed only slightly worse than the field camera measurement. 

- `get_test_wave.m` Function to create arbitrary Pulseq block from all the different test waves used in the GIRF acquisition.




#### Data Files

- `girf_skope2.seq` contains the exact Pulseq sequences used to measure the GIRF in this work.

- `girf3d2_XXXX.seq` contains the Pulseq sequences using the thinslice method, where XXXX are different possble test waveforms.

- `girf_skope_waves.mat` This contains all the test gradient waveforms used for calculating the GIRF

- `bipolars.mat` and `wave_save_v2.mat` contian pre-computed test waves for the GIRF acquisition.



### Installation

Matlab scripts needs Pulseq to be in your path.

Python script modules should all be available on common package managers.