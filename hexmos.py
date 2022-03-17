# Define a few functions

def create_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-ms", "--mos-size", type=float, default=5.,
                   help="Size of the square region within which to distribute the pointings (deg). "
                        " Default is 5 deg.")
    p.add_argument("-bw", "--beam-width", type=float, default=1.,
                   help="Primary beam width (FWHM in deg; a Gaussian shape is assumed). "
                        "Default is 1 deg.")
    p.add_argument("-gs", "--grid-spacing", type=float, default=0.58,
                   help="Spacing of the hexagonal mosaic grid (deg). "
                        "Default is 0.58 deg = 1 deg / sqrt(3) = Nyquist sampling for default beam.")
    return p

def beam(X,Y,FWHM):
  # Gaussian beam
  return np.exp(-4*np.log(2)*((X)**2+(Y)**2)/FWHM**2)


# Import modules

import numpy as np
from matplotlib import pyplot as plt
import argparse,sys


# Read settings from command line

args = create_parser().parse_args([a for a in sys.argv[1:]])
ms = args.mos_size
bw = args.beam_width
gs = args.grid_spacing


# Print settings

print('#')
print('# Will make a mosaic with the following specs:')
print('    - mosaic size = {0:1.2f} deg'.format(ms))
print('      (size of the square region within which to distribute the pointings)')
print('    - primary beam FWHM = {0:1.2f} deg'.format(bw))
print('      (assuming a Gaussian beam)')
print('    - hexagonal mosaic grid spacing = {0:1.2f} deg = FWHM / {1:1.2f}'.format(gs, bw/gs))
print('#')


# Hardcoded settings

pix = bw/100 # pixel size of mosaic map
show_point = True


# Go!

coords=np.atleast_2d(np.arange(-1.2*ms/2,1.2*ms/2+pix,pix))
x=coords.repeat(coords.shape[1],axis=0)
coords=np.transpose(coords)
y=coords.repeat(coords.shape[0],axis=1)
field=np.zeros(x.shape)

fig=plt.figure(figsize=(6,6))
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.98,top=0.98,wspace=0.,hspace=0.)
ax=fig.add_subplot(111, aspect='equal')

yp=np.arange(0,ms/2,gs)
xp=np.arange(0,ms/2,gs*np.sqrt(3))
Np=0
for yy in yp:
  for xx in xp:
    if not xx and not yy:
      if show_point: plt.plot([0],[0],'k.')
      field=field+beam(x,y,bw)**2
      Np+=1
    elif not xx and yy:
      if show_point: plt.plot([0],[-yy],'k.')
      field=field+beam(x,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([0],[+yy],'k.')
      field=field+beam(x,y+yy,bw)**2
      Np+=1
    elif xx and not yy:
      if show_point: plt.plot([-xx],[0],'k.')
      field=field+beam(x-xx,y,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[0],'k.')
      field=field+beam(x+xx,y,bw)**2
      Np+=1
    else:
      if show_point: plt.plot([-xx],[-yy],'k.')
      field=field+beam(x-xx,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([-xx],[+yy],'k.')
      field=field+beam(x-xx,y+yy,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[-yy],'k.')
      field=field+beam(x+xx,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[+yy],'k.')
      field=field+beam(x+xx,y+yy,bw)**2
      Np+=1

yp=yp+gs/2
yp=yp[yp<=ms/2]
xp=xp+gs*np.sqrt(3)/2
xp=xp[xp<=ms/2]
for yy in yp:
  for xx in xp:
    if not xx and not yy:
      if show_point: plt.plot([0],[0],'k.')
      field=field+beam(x,y,bw)**2
      Np+=1
    elif not xx and yy:
      if show_point: plt.plot([0],[-yy],'k.')
      field=field+beam(x,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([0],[+yy],'k.')
      field=field+beam(x,y+yy,bw)**2
      Np+=1
    elif xx and not yy:
      if show_point: plt.plot([-xx],[0],'k.')
      field=field+beam(x-xx,y,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[0],'k.')
      field=field+beam(x+xx,y,bw)**2
      Np+=1
    else:
      if show_point: plt.plot([-xx],[-yy],'k.')
      field=field+beam(x-xx,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([-xx],[+yy],'k.')
      field=field+beam(x-xx,y+yy,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[-yy],'k.')
      field=field+beam(x+xx,y-yy,bw)**2
      Np+=1
      if show_point: plt.plot([+xx],[+yy],'k.')
      field=field+beam(x+xx,y+yy,bw)**2
      Np+=1

field=1./np.sqrt(field)
fmin=field.min()


# Print results

print('# Result:')
print('    {0:d} pointings, min noise = {1:.2f} x single-pointing noise'.format(Np,fmin))
print('    mosaic area = {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<1.1*fmin).sum())*pix**2,1.1))
print('                  {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<np.sqrt(2)*fmin).sum())*pix**2,np.sqrt(2)))
print('                  {0:.1f} deg^2 below {1:.2f}x min noise'.format(float((field<2*fmin).sum())*pix**2,2))
print('#')


# Plot

ax.contour(x,y,field, [1.1*fmin, np.sqrt(2)*fmin, 2*fmin], colors='k')
plt.xlabel('X (deg)')
plt.ylabel('Y (deg)')
plt.show()
