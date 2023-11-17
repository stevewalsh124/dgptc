## dgptc
Deep Gaussian processes (DGPs) to quantify uncertainty of computer model output related to cosmological redshift (`1D_examples`) tropical cyclone precipitation forecasts (`2D_` and `3D_examples`).

### `1D_examples`

`1D_real_study_full_emuspace.R` is the main script to obtain uncertainty quantification (UQ) using DGPs for the cosmological dark matter power spectrum, given a particular batch of computer model runs (there are 111 batches of model runs). These results are calculated on the emulation space, $\mathcal{P}(k) = \log_{10}(k^{1.5} P(k) / (2 \pi^2))$ for the full range of $k$ values. In order to get the results for each of the 111 batches, you can use the corresponding shell script and run the following in the command line (you will have to update your own ARC info within the shell script):

`for i in $(seq 1 111); do sbatch 1D_real_study_full_emuspace.sh $i; done`.

For more details on the background of this data as well as the subsequent modeling done, see the third chapter of my [dissertation](https://vtechworks.lib.vt.edu/bitstream/handle/10919/115494/Walsh_SA_D_2023.pdf?sequence=1&isAllowed=y). After getting the results from each batch of runs, you can compare how the posterior mean compares to that of the cosmic emulator ([Cosmic Emu](https://github.com/lanl/CosmicEmu#readme)) using the `all_UQ_fullemu.R` script. In order to use the Cosmic Emu, you will need an `xstar` file; my `xstar.R` creates a corresponding file based on a particular input design of interest. 
