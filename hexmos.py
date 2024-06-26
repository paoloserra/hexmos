# Define a few functions
def create_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-ms", "--mos-size", type=str, default='5,5',
                   help="X,Y size of the square region within which to distribute the pointings (deg)."
                        " Default is 5,5 deg.")
    p.add_argument("-bw", "--beam-width", type=float, default=1.,
                   help="Primary beam width (FWHM in deg; a Gaussian shape is assumed). "
                        "Default is 1 deg.")
    p.add_argument("-gs", "--grid-spacing", type=float, default=0.58,
                   help="Spacing of the hexagonal mosaic grid (deg). "
                        "Default is 0.58 deg = 1 deg / sqrt(3) = Nyquist sampling for default beam.")
    p.add_argument("-ra", "--right-ascension", type=str, default='12:0:0.0',
                   help="Right ascension of the mosaic centre in hh:mm:ss.s. Default is 12:0:0.0.")
    p.add_argument("-dec", "--declination", type=str, default='30:0:0.0',
                   help="Declination of the mosaic centre in dd:mm:ss.s. Default is 30:0:0.")
    p.add_argument("-rot", "--rotate", type=float, default=0.0,
                   help="Rotate pointing pattern about the mosaic centre by this angle (north through east)."
                   " Default is 0 deg.")
    p.add_argument("-sp", "--show-pointings", action="store_true",)
    p.add_argument("-sc", "--save-coords", type=str,
                   help="Save the coordinates of the pointings to an ASCII file. Give the file name here."
                   " The default is to not save the coordinates to a file.")
    p.add_argument("-skip", "--skip-pointings", type=int, nargs='+', required = False,
                   help="Space separated list of pointings ID to be skipped when building the mosaic."
                   " Default is to include all poitings.")
    p.add_argument("-aup", "--add-user-plot", action="store_true",
                   help="Include additional elements in the plot through the python module called"
                   " add_user_plot.py located in the CWD. The module must define the single function"
                   " add_to_plot(axis,wcs), where 'axis' is a matplotlib.axes.Axes object and"
                   " 'wcs' is a astropy.wcs.WCS object.")
    return p

def beam(X,Y,FWHM):
  # Gaussian beam
  return np.exp(-4*np.log(2)*((X)**2+(Y)**2)/FWHM**2)


# Import modules
import numpy as np
from matplotlib import pyplot as plt
import argparse,sys
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.wcs as WCS
from astropy.io import fits


# Fix ArgumentParser bug with negative arguments
if '-dec=' in sys.argv[1:]:
    decsign = ''
elif '-dec' in sys.argv[1:]:
    decvalind = sys.argv.index('-dec')+1
    if sys.argv[decvalind][0] == '-':
        sys.argv[decvalind] = sys.argv[decvalind][1:]
        decsign = '-'
    else: decsign = ''
else: decsign = ''


# Read settings from command line
args = create_parser().parse_args([a for a in sys.argv[1:]])
ms = list(map(float,args.mos_size.split(',')))
bw = args.beam_width
gs = args.grid_spacing
sp = args.show_pointings
sc = args.save_coords
rot = args.rotate
skip = args.skip_pointings
aup = args.add_user_plot

mos_cen_ra = args.right_ascension
mos_cen_dec = decsign+args.declination
if not skip:
    skip = []


# Print settings
print('#')
print('# Will make a mosaic with the following specs:')
print('#    - mosaic centre RA(h:m:s), Dec(d:m:s) = {0:s} , {1:s}'.format(mos_cen_ra, mos_cen_dec))
print('#    - mosaic size = {0:1.2f} deg x {1:1.2f} deg (RA axis x Dec axis)'.format(ms[0], ms[1]))
print('#      (size of the region within which to distribute the pointings)')
print('#    - mosaic rotation = {0:.1f} deg north through east'.format(rot))
print('#    - primary beam FWHM = {0:1.2f} deg'.format(bw))
print('#      (assuming a Gaussian beam)')
print('#    - hexagonal mosaic grid spacing = {0:1.2f} deg = FWHM / {1:1.2f}'.format(gs, bw/gs))
print('#')


# Hardcoded settings
pix = bw/100 # pixel size of mosaic map


# Go!


# Create the arrays with the X and Y coordinates
rot *= -np.pi / 180
array_size = 1.2*np.max((
                np.abs(np.cos(rot) * ms[0]) + np.abs(np.sin(rot) * ms[1]),
                np.abs(np.sin(rot) * ms[0]) + np.abs(np.cos(rot) * ms[1])))
x = np.atleast_2d(np.arange(-array_size/2, array_size/2+pix, pix))
y = np.atleast_2d(np.arange(-array_size/2, array_size/2+pix, pix))
y = np.transpose(y)
x = x.repeat(y.shape[0],axis=0)
y = y.repeat(x.shape[1],axis=1)


# Create sky image with axis0=Y, axis1=X
field=np.zeros((y.shape[0],x.shape[1]))


# Creat list of pointings on hexagonal grid
points = []
yp=np.arange(0,ms[1]/2,gs)
xp=np.arange(0,ms[0]/2,gs*np.sqrt(3))
for yy in yp:
  for xx in xp:
    if not xx and not yy:
      if not [0, 0] in points: points.append([0, 0])
    elif not xx and yy:
      if not [0, -yy] in points: points.append([0, -yy])
      if not [0, +yy] in points: points.append([0, +yy])
    elif xx and not yy:
      if not [-xx, 0] in points: points.append([-xx, 0])
      if not [+xx, 0] in points: points.append([+xx, 0])
    else:
      if not [-xx, -yy] in points: points.append([-xx, -yy])
      if not [-xx, +yy] in points: points.append([-xx, +yy])
      if not [+xx, -yy] in points: points.append([+xx, -yy])
      if not [+xx, +yy] in points: points.append([+xx, +yy])
yp=yp+gs/2
yp=yp[yp<=ms[1]/2]
xp=xp+gs*np.sqrt(3)/2
xp=xp[xp<=ms[0]/2]
for yy in yp:
  for xx in xp:
    if not xx and not yy:
      if not [0, 0] in points: points.append([0, 0])
    elif not xx and yy:
      if not [0, -yy] in points: points.append([0, -yy])
      if not [0, +yy] in points: points.append([0, +yy])
    elif xx and not yy:
      if not [-xx, 0] in points: points.append([-xx, 0])
      if not [+xx, 0] in points: points.append([+xx, 0])
    else:
      if not [-xx, -yy] in points: points.append([-xx, -yy])
      if not [-xx, +yy] in points: points.append([-xx, +yy])
      if not [+xx, -yy] in points: points.append([+xx, -yy])
      if not [+xx, +yy] in points: points.append([+xx, +yy])
points = np.array(points)

# Rotate pointings coordinates
points = np.concatenate((
                         np.cos(rot) * points[:,0:1] - np.sin(rot) * points[:,1:2],
                         np.sin(rot) * points[:,0:1] + np.cos(rot) * points[:,1:2]),
                         axis = 1)

# Add beams to sky image
mos_cen = SkyCoord('{0:s} {1:s}'.format(mos_cen_ra,mos_cen_dec), frame='icrs', unit=(u.hourangle, u.deg))
if sc:
    f = open(sc, 'w')
    f.write('# {0:>5s} {1:>8s} {2:>8s} {3:>8s} {4:>8s}\n'.format('ID', 'Xoffset', 'Yoffset', 'RA', 'Dec'))
    f.write('# {0:>5s} {1:>8s} {2:>8s} {3:>8s} {4:>8s}\n'.format('', '(deg)', '(deg)', '(deg)', '(deg)'))
Np=0
radecs = []
for pp in points:
    [xx, yy] = pp
    radec = mos_cen.spherical_offsets_by(xx*u.deg, yy*u.deg)
    if Np in skip:
        if sc: f.write('#  {0:5d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f}\n'.format(Np, xx, yy, radec.ra.value, radec.dec.value))
    else:
        field=field+beam(x-xx, y-yy, bw)**2
        radecs.append([Np,radec])
        if sc: f.write('  {0:5d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f}\n'.format(Np, xx, yy, radec.ra.value, radec.dec.value))
    Np+=1
if sc:
    f.close()


# Convert sky image into noise pattern
field = 1./np.sqrt(field)
field = field[:,::-1] # to account for RA axis direction
fmin = field.min()


# Print results
print('# Result:')
print('#    {0:4d} pointings, min noise = {1:.2f} x single-pointing noise'.format(len(radecs),fmin))
print('#                    max effective integration = {0:.2f} x single-pointing integration'.format(1./fmin**2))
print('#    mosaic area = {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<1.1*fmin).sum())*pix**2,1.1))
print('#                  {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<np.sqrt(2)*fmin).sum())*pix**2,np.sqrt(2)))
print('#                  {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<2*fmin).sum())*pix**2,2))
print('#')

field[field > 100 * np.nanmin(field)] = np.nan

# Create FITS file for plotting and save it to disk
hdu = fits.PrimaryHDU(field)
hdu.header['naxis']   = 2
hdu.header['naxis1']  = field.shape[1]
hdu.header['naxis2']  = field.shape[0]
hdu.header['crpix1']  = (field.shape[1]+1)/2
hdu.header['crpix2']  = (field.shape[0]+1)/2
hdu.header['crval1']  = mos_cen.ra.deg
hdu.header['crval2']  = mos_cen.dec.deg
hdu.header['cdelt1']  = -pix
hdu.header['cdelt2']  = +pix
hdu.header['ctype1']  = 'RA---TAN'
hdu.header['ctype2']  = 'DEC--TAN'
hdu.header['equinox'] = 2000
hdu.writeto('snake.fits', overwrite=True)


# Initialise figure
fig = plt.figure(figsize=(7,7))
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.98,top=0.98,wspace=0.,hspace=0.)
ax = fig.add_subplot(111, aspect='equal', projection=WCS.WCS(hdu.header))


# Add pointings and sensitivity contours to figure
if sp:
    for pp in radecs:
        x,y = WCS.utils.skycoord_to_pixel(pp[1], WCS.WCS(hdu.header))
        ax.plot(x, y, marker='.', ms=0.5, color='0.5')
        ax.text(x, y, pp[0], ha='center', va='center', fontsize=4, zorder=100)
ax.contour(field, [1.1*fmin, np.sqrt(2)*fmin, 2*fmin], colors='0.5', linewidths=[0.5,], zorder=10)

if aup:
    import add_user_plot as aup
    aup.add_to_plot(ax, WCS.WCS(hdu.header))

# Finalise figure
plt.xlabel('RA (J2000)')
plt.ylabel('Dec (J2000)')
plt.grid()
plt.show()
