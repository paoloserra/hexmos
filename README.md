# hexmos
Interferometric mosaic planning tool.


It distributes Gaussian primary beams of a given FWHM (`-bw`) on a hexagonal mosaic grid of a given size (`-ms`) and spacing (`-gs`), and returns useful info on the resulting mosaic.

```
usage: hexmos.py [-h] [-ms MOS_SIZE] [-bw BEAM_WIDTH] [-gs GRID_SPACING]

optional arguments:
  -h, --help            show this help message and exit
  -ms MOS_SIZE, --mos-size MOS_SIZE
                        Size of the square region within which to distribute the pointings (deg). Default is 5 deg.
  -bw BEAM_WIDTH, --beam-width BEAM_WIDTH
                        Primary beam width (FWHM in deg; a Gaussian shape is assumed). Default is 1 deg.
  -gs GRID_SPACING, --grid-spacing GRID_SPACING
                        Spacing of the hexagonal mosaic grid (deg). Default is 0.58 deg = 1 deg / sqrt(3) = Nyquist sampling for default beam.
```
