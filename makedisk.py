"""Module with methods for use in constructing a T Tauri disk cube for the FISICA 
simulator.

The methods are:
"""

from __future__ import absolute_import

import math
import numpy as np
import scipy.interpolate as interpolate
import scipy.constants as sc
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt

import common.commonobjects as co
import makecube
import writefits

def makeDisk(wn_min, wn_max, wn_step, beta=1.5):

    epsilon = 1.0
    rtap = 50.0
    h0 = 10.0
    r0 = 100.0
    disk_beta = 1.1
    rinner = 0.07

#    ns = 2000
    ns = 400
#    max_r = 100.0
    rmax = 200.0

    x = np.linspace(-rmax, rmax, ns)
    y = np.linspace(-rmax, rmax, ns)
    z = np.linspace(-rmax, rmax, ns)
    
    xm, ym, zm = np.meshgrid(x, y, z)

    theta = np.deg2rad(40.0)
    # xd, yd, zd are coords in disk frame disk centre at 0,0,ns/2, disk
    # pole along z axis
    xd = xm
    yd = ym * np.cos(theta) + zm * np.sin(theta)
    zd = -ym * np.sin(theta) + zm * np.cos(theta)

    r2d = np.sqrt(xd**2 + yd**2)
    r3d = np.sqrt(xd**2 + yd**2 + zd**2)

    column_density = r2d**(-epsilon) * np.exp(-(r2d / rtap)**(2 - epsilon))
    column_density[r2d < rinner] = rinner**(-epsilon) * np.exp(-(rinner / rtap)**(2 - epsilon))

    scale_height = h0 * (r2d / r0)**disk_beta 
    scale_height[r2d < rinner] = h0 * (rinner / r0)**disk_beta 

    print r2d.shape
    print r2d[:,:,0]

    # relative density = n_midplane * exp (-h/scale_h)
    density = (column_density / (2 * scale_height)) * \
      np.exp(-(np.abs(zd) / scale_height))
    density[r3d > rmax - 10] = 0.0 

    # disk mass, assume all gas, kilograms
    solar_mass = 1.988e30
    mdisk = 3e-2 * solar_mass
    
    # following assumes disk all H, 20% inaccurate
    mh = 1.67e-27

    # total mass of uncorrected disk model
    disk_total = np.sum(density)

    mass_per_voxel = mdisk / disk_total
    au = 1.496e11
    voxel_side = np.abs(x[1] - x[0]) * au
    voxel_vol = voxel_side**3
    nh_per_voxel = mass_per_voxel / (voxel_vol * mh) 
    
    # convert disk model to nh per cubic metre
    density *= nh_per_voxel

    check_mass = np.sum(density) * voxel_vol * mh
    print 'check mass', check_mass

    # now convert to Av (tau_v) assuming Av~1 at N_H~1e21cm-2 = 1e25m-2
    a_v = density * voxel_side / 1.0e25

    # now loop through wavenumbers of interest
    wn_list = np.linspace(wn_min, wn_max,
      int(np.rint((wn_max-wn_min)/wn_step)) + 1)
    bandwidth = wn_step * 3.0e10

    # temperature (profile numbers are log10 r, log10 T)
    temp_profile = [(-1.5, 2.85),
      (-0.92, 2.477),
      (4.3e-3, 2.0),
      (0.602, 1.602),
      (1.3, 1.3),
      (1.6, 1),
      (2.0, 1),
      (3.0, 0.8)]
    r = []
    t = []
    for point in temp_profile:
        r.append(point[0])
        t.append(point[1])
    f = interpolate.interp1d(r, t, 'linear')
    print f([-1, 0, 1, 2, 2.9])
    temperature = 10**f(np.log10(r2d))

    cubeshape = wn_list.shape + (ns, ns,)
    emission = np.zeros(cubeshape, np.float)

    for iwn, wn in enumerate(wn_list):
#    for iwn, wn in enumerate([150.0]):
        print wn

        bb = makecube.black_body(temperature, wn)

        # V~500nm
        tau = a_v * (500e-9 * wn * 100.0)**beta

        # tau_integral measured from front of cube
        tau_integral = np.cumsum(tau, axis=2)

        # emission in Jy
        emission[iwn, :, :] = 1e26 * np.sum(
          bb * (1.0 - np.exp(-tau)) * np.exp(-tau_integral), axis=2)

    # write the cube to FITS
#    axis1 = co.Axis(data=x, title='x', units='arcsec')
#    axis2 = co.Axis(data=y, title='y', units='arcsec')
#    axis3 = co.Axis(data=z, title='z', units='arcsec')

    axis0 = co.Axis(data=wn_list, title='frequency', units='cm-1')
    axis1 = co.Axis(data=y, title='x', units='AU')
    axis2 = co.Axis(data=z, title='y', units='AU')

    cube = co.Image(data=emission, axes=[axis0, axis1, axis2])
    cubewriter = writefits.WriteFITSCube(cube, 'makedisk.fits')
    cubewriter.run()

def makeDiskImage(max_baseline, m1_diameter=2.0, disk_cube='makedisk.fits',
  distance=100.0):
    """Method to construct a 2d image with spatial characteristics 
    suitable for the FISICA simulator from an input 2d FITS image.

    Keyword arguments:
    max_baseline  -- maximum baseline of the simulation [m]
    wn_min        -- lowest frequency of the simulation [cm-1]
    wn_max        -- highest frequency of the simulation [cm-1]
    m1_diameter   -- diameter of collector primary mirror [m] (default 2)

    Returns:
    A commonobjects.Image object containing the input image mapped onto
    a grid covering the primary beam (to first 0) of the interferometer
    and sampled at the Nyquist frequency of the longest baseline / highest
    frequency
    """

    hdulist = pyfits.open(disk_cube)
    print hdulist.info()
    prihdr = hdulist[0].header
    print repr(prihdr)
    fits_data = hdulist[0].data

    if np.any(np.isnan(fits_data)):
        print 'image array contains NaN values, setting these to 0'
        fits_data[np.isnan(fits_data)] = 0.0

    # Note that, like C (and unlike FORTRAN), Python is 0-indexed and 
    # the indices have the slowest axis first and fastest changing axis
    # last; i.e. for a 2-D image, the fast axis (X-axis) which 
    # corresponds to the FITS NAXIS1 keyword, is the second index.

    datashape = np.shape(fits_data)
    print datashape
    nwn = datashape[0]
    ny = datashape[1]
    nx = datashape[2]

    try:
        object_name = prihdr['OBJECT']
    except:
        object_name = 'unknown'

    cdelt1 = prihdr['CDELT1']
    cdelt2 = prihdr['CDELT2']
    cdelt3 = prihdr['CDELT3']

    crval1 = prihdr['CRVAL1']
    crval2 = prihdr['CRVAL2']
    crval3 = prihdr['CRVAL3']

    crpix1 = prihdr['CRPIX1']
    crpix2 = prihdr['CRPIX2']
    crpix3 = prihdr['CRPIX3']

    ctype1 = prihdr['CTYPE1']
    ctype2 = prihdr['CTYPE2']
    ctype3 = prihdr['CTYPE3']

    cunit1 = prihdr['CUNIT1']
    assert 'AU' in cunit1 
    cunit2 = prihdr['CUNIT2']
    assert 'AU' in cunit2 
    cunit3 = prihdr['CUNIT3']
    assert 'HZ' in cunit3

    bunit = prihdr['BUNIT']

    # + 1 is in there because FITS arrays are 1-based
    fits_axis1 = crval1 + (np.arange(nx) + 1 - crpix1) * cdelt1
    fits_axis2 = crval2 + (np.arange(ny) + 1 - crpix2) * cdelt2
    fits_axis3 = crval3 + (np.arange(nwn) + 1 - crpix3) * cdelt3

    assert 'JY/PIXEL' in bunit.upper()

    hdulist.close()

    # calculate pixel size of desired cube, should be as small as
    # the Nyquist freq of the longest baseline / highest freq
    wn_max = np.max(fits_axis3) / 3e10
    print 'wn_max', wn_max
    target_pixsize = 1.0 / (wn_max * 100.0 * max_baseline)
    target_pixsize /= 2.0
    target_pixarea = target_pixsize**2
    target_pixsize = np.rad2deg(target_pixsize) * 3600
    print 'target pixel size', target_pixsize, 'arcsec'

    # now get size of target image to be created. This will be 
    # square, big enough to cover the primary beam, and centred at 
    # the FITS data array centre
    wn_min = np.min(fits_axis3) / 3e10
    print 'wn_min', wn_min
    primary_beam_r = 1.22 / (wn_min * 100 * m1_diameter)
    primary_beam_r = np.rad2deg(primary_beam_r) * 3600 
    target_extent = 2 * primary_beam_r
    print 'primary beam diameter', primary_beam_r, 'arcsec'

    tnpix = int(target_extent / target_pixsize)

    trpix = int(tnpix / 2)
    trval = 0
    tdelt = target_pixsize
    target_axis = trval + (np.arange(tnpix, dtype=np.float) - trpix) * \
      tdelt 

    print 'target image will be', tnpix, 'x', tnpix, 'square'

    # convert fits spatial axes to arcsec
    # 1 AU at 1 parsec subtends 1 arcsec
    apparent_axis1 = fits_axis1 / distance
    apparent_axis2 = fits_axis2 / distance
    # pixel area in sr
    fits_pixarea = np.abs(
      np.deg2rad((apparent_axis1[1] - apparent_axis1[0]) / 3600) * \
      np.deg2rad((apparent_axis2[1] - apparent_axis2[0]) / 3600))

    smooth = abs(target_pixsize / (apparent_axis1[1] - apparent_axis1[0]))
    print 'smoothing parameter', smooth

    target_data = np.zeros([nwn, tnpix, tnpix], dtype=np.float)
    for iwn, wn in enumerate(fits_axis3):

        # interpolate from FITS image onto target image
        # ..RectBivariateSpline requires axis1 and axis2 be monotonically
        # ..increasing, flip if need to
#    axis1_mult = np.sign(magnified_axis1[1] - magnified_axis1[0])
#    magnified_axis1 *= axis1_mult

#    axis2_mult = np.sign(magnified_axis2[1] - magnified_axis2[0])
#    magnified_axis2 *= axis2_mult

        print np.shape(apparent_axis1), np.shape(apparent_axis2), np.shape(fits_data)
        print '..setting up the interpolation'
        interp = interpolate.RectBivariateSpline(apparent_axis2, 
          apparent_axis1, fits_data[iwn,:,:], s=0)
        print '..interpolating'
        target_data[iwn,:,:] = interp(target_axis, target_axis, grid=True)

        print np.max(fits_data[iwn,:,:]), np.max(target_data[iwn,:,:])

    # reject obviously rubbish spline interpolations
    target_data[target_data < 0] = 0.0

    # multiply interpolated image by ratio of pixel areas to make
    # units brightness/pixel for the new pixel size
    print target_pixarea / fits_pixarea

    target_data *= (target_pixarea / fits_pixarea)

    axis0 = co.Axis(data=fits_axis3, title='frequency', units='Hz')
#    axis1 = co.Axis(data=fits_axis1, title='RA', units='arcsec') 
#    axis2 = co.Axis(data=fits_axis2, title='Dec', units='arcsec') 
#    cube = co.Image(fits_data, title=object_name, axes=[axis0, axis1, axis2],
#      units='Jy/pixel')
    axis1 = co.Axis(data=target_axis, title='RA', units='arcsec') 
    axis2 = co.Axis(data=target_axis, title='Dec', units='arcsec') 
    cube = co.Image(target_data, title=object_name, axes=[axis0, axis1, axis2],
      units='Jy/pixel')

    cubewriter = writefits.WriteFITSCube(cube, 'makediskimage.fits')
    cubewriter.run()
