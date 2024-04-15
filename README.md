# hexmos
Interferometric mosaic planning tool.


It distributes Gaussian primary beams of a given FWHM (`-bw`) on a hexagonal mosaic grid of a given size (`-ms`) and spacing (`-gs`), and returns useful info on the resulting mosaic.

Users can specify the centra RA (`-ra`) and Dec (`-dec`), and a rotation angle about the centre (`-rot`).

The RA,Dec of the pointings can be saved to an ASCII file (`-sc`).

It is possible to modify a mosaic set up exclude some pointgs (`-skip`).

```
usage: hexmos.py [-h] [-ms MOS_SIZE] [-bw BEAM_WIDTH] [-gs GRID_SPACING] [-ra RIGHT_ASCENSION] [-dec DECLINATION] [-rot ROTATE] [-sc SAVE_COORDS] [-skip SKIP_POINTINGS [SKIP_POINTINGS ...]]

optional arguments:
  -h, --help            show this help message and exit
  -ms MOS_SIZE, --mos-size MOS_SIZE
                        X,Y size of the square region within which to distribute the pointings (deg). Default is 5,5 deg.
  -bw BEAM_WIDTH, --beam-width BEAM_WIDTH
                        Primary beam width (FWHM in deg; a Gaussian shape is assumed). Default is 1 deg.
  -gs GRID_SPACING, --grid-spacing GRID_SPACING
                        Spacing of the hexagonal mosaic grid (deg). Default is 0.58 deg = 1 deg / sqrt(3) = Nyquist sampling for default beam.
  -ra RIGHT_ASCENSION, --right-ascension RIGHT_ASCENSION
                        Right ascension of the mosaic centre in hh:mm:ss.s. Default is 12:0:0.0.
  -dec DECLINATION, --declination DECLINATION
                        Declination of the mosaic centre in dd:mm:ss.s. Default is 30:0:0.
  -rot ROTATE, --rotate ROTATE
                        Rotate pointing pattern about the mosaic centre by this angle (north through east). Default is 0 deg.
  -sc SAVE_COORDS, --save-coords SAVE_COORDS
                        Save the coordinates of the pointings to an ASCII file. Give the file name here. The default is to not save the coordinates to a file.
  -skip SKIP_POINTINGS [SKIP_POINTINGS ...], --skip-pointings SKIP_POINTINGS [SKIP_POINTINGS ...]
                        Space separated list of pointings ID to be skipped when building the mosaic. Default is to include all poitings.
```
