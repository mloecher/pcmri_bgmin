### Installation

Matlab scripts needs Pulseq to be in your path.

To run the python code, you must install Gropt.  Currently the optimization routine is in a test pip repository and can be installed by running:

`pip install -U --extra-index-url https://test.pypi.org/simple/ gropt2`

If you get an error with including the fftw library, you may need to add: `os.add_dll_directory('fftw-path')` or installing the `pyfftw` package also seems to work.



### Files

- `simple_waves.ipynb` is a Python notebook that was used to generate all the waveforms used in the sliding TE experiments.  It generates waveforms with different echo times (minimization window positions), as well as different slice orientations (by rotating the GIRF being used).  

- `loop_seqs_v3.m` Generates Pulseq files used to conduct the sliding TE exepriments using the waveforms generated in `simple_waves.ipynb`.

- `single_seq_r2.m` Generates Pulseq files used to measure flow in vivo, using the GrOpt generated waveforms as arbitrarily shaped gradient blocks.  Requires recording of ECG and respiratory signals for retrospective binning.