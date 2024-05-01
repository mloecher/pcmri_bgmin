# pcmri_bgmin
Code to minimize PC-MRI background phase based on GIRF predictions.

## Details
This code gives demonstrations in two parts:

- GIRF/

Contains acquisition and reconstruction code to generate the GIRFs used in this work, including Pulseq files that can be run on a scanner to measure your own GIRFs

- Waveforms/

Contains example Python code to use GrOpt to generate PC-MRI background minimizing waveforms from a set of parameters and a measured GIRF.