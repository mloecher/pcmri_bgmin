### Files

- `girf_utils.py` contains example functions to generate a GIRF from the measurements acquired from the field camera measurements.
- `girf_v2.m` contains code to generate a GIRF measurement sequence for use with a Skope field camera.
- `girf_skope2.seq` contains the exact Pulseq sequenced used to measure the GIRF in this work.
- `girf3d_2.m` is a Matlab script for generating a Pulseq sequence to measure a 3D GIRF using the thin slice method.  This was not thoroughly tested (because we switched to the field camera), however initial tests showed that it performed only slightly worse than the field camera measurement. 