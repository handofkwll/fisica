from __future__ import absolute_import

import numpy as np
import scipy.interpolate as interpolate
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt

import common.commonobjects as co
from skygenerator import BB_spectrum
import writefits


class AddCubes(object):

    def __init__(self, in_cubes=[], cubename='total_cube.fits'):
        self.in_cubes = in_cubes
        self.cubename = cubename

    def run(self):
        print 'AddCubes.run'

        print 'template FITS file', self.in_cubes[0]
        hdulist = pyfits.open(self.in_cubes[0])
        skydata = hdulist[0].data
#        hdulist.close()

        for file in self.in_cubes[1:]:
            print 'adding file', file
            new_hdulist = pyfits.open(file)
            new_data = new_hdulist[0].data
            print new_data.shape, skydata.shape
            skydata += new_data
            new_hdulist.close()

        # write the cube to FITS
        print 'writing', self.cubename
        hdulist[0].data = skydata
        hdulist.writeto(self.cubename, clobber=True)
        hdulist.close()


class MakeCube(object):
    """Class to make a cube from a 2-d image crossed with a spectrum. 
    """

    def __init__(self, image, spectrum, cubename='cube.fits'):
        self.image = image
        self.spectrum = spectrum
        self.cubename = cubename

    def run(self):
        print 'MakeCube.run'

        # the cube shape is just the image replicated through the spectral
        # dimension - [wn, dec, ra]
        cubeshape = self.spectrum.data.shape + self.image.data.shape
        cube = np.zeros(cubeshape, np.float)

        # normalise the image plane so that integrated flux from the 
        # cube is the spectrum. If the spectrum is in W/m2/Hz then cube
        # units will be W/m2/Hz/voxel.
        self.image.data /= np.sum(self.image.data)
      
        # populate each plane of the cube in turn
        for iwn,wn in enumerate(self.spectrum.axis.data):
            cube[iwn,:,:] = self.image.data * self.spectrum.data[iwn]

        # write the cube to FITS
        cube = co.Image(data=cube, axes=[self.spectrum.axis] + self.image.axes)
        cubewriter = writefits.WriteFITSCube(cube, self.cubename)
        cubewriter.run() 
        

class MakeImage(object):
    """
    """

    def __init__(self, fits_image, max_baseline, wn_nominal, wn_min, wn_max,
      magnification, m1_diameter=2.0, plotname=None, verbose=False):
        assert wn_max > wn_min

        self.fits_image = fits_image
        self.max_baseline = max_baseline
        self.wn_nominal = wn_nominal
        self.wn_max = wn_max
        self.wn_min = wn_min
        self.magnification = magnification
        self.m1_diameter = m1_diameter
        self.plotname = plotname
        self.verbose = verbose

    def run(self):
        print 'MakeImage.run'

        hdulist = pyfits.open(self.fits_image)
        print hdulist.info()
        prihdr = hdulist[0].header
        if self.verbose:
            print repr(prihdr)
        fits_data = hdulist[0].data

        if np.any(np.isnan(fits_data)):
            print 'image array contains NaN values, setting these to 0'
            fits_data[np.isnan(fits_data)] = 0.0

        # Note that, like C (and unlike FORTRAN), Python is 0-indexed and 
        # the indices have the slowest axis first and fastest changing axis
        # last; i.e. for a 2-D image, the fast axis (X-axis) which 
        # corresponds to the FITS NAXIS1 keyword, is the second index.

        # get rid of degenerate axes (freq), verify data are 2d and swap
        # axes - fits_data ends with shape [AXIS2, AXIS1]
        fits_data = np.squeeze(fits_data)
        assert len(fits_data.shape) == 2
        # for Fomalhaut don't swap axes - looks like Fomalhaut data is in C order
        fits_data = np.swapaxes(fits_data, 0, 1)

        datashape = np.shape(fits_data)
        nx = datashape[1]
        ny = datashape[0]

        try:
            object_name = prihdr['OBJECT']
        except:
            object_name = 'unknown'

        # pixel increments in degrees
        cdelt1 = prihdr['CDELT1']
        cdelt2 = prihdr['CDELT2']

        crval1 = prihdr['CRVAL1']
        crval2 = prihdr['CRVAL2']

        crpix1 = prihdr['CRPIX1']
        crpix2 = prihdr['CRPIX2']

        ctype1 = prihdr['CTYPE1']
        ctype2 = prihdr['CTYPE2']

        # spatial coords are not stored in Fomalhaut FITS file
        # but look like they may be degrees
        try:
            cunit1 = prihdr['CUNIT1']
        except:
            print "CUNIT1 not set, assuming 'DEG'"
            cunit1 = 'DEG'

        try:
            print "CUNIT2 not set, assuming 'DEG'"
            cunit2 = prihdr['CUNIT2']
        except:
            cunit2 = 'DEG'

        bunit = prihdr['BUNIT']

        # calculate axes of FITS image in arcsec
        if 'DEG' in cunit1.upper():
            cdelt1 *= 3600.0
        else:
            raise Exception, 'cannot handle CUNIT1=%s' % cunit1

        if 'DEG' in cunit2.upper():
            cdelt2 *= 3600.0
        else:
            raise Exception, 'cannot handle CUNIT2=%s' % cunit2

        # + 1 is in there because FITS arrays are 1-based
        fits_axis1 = (np.arange(nx) + 1 - crpix1) * cdelt1
        fits_axis2 = (np.arange(ny) + 1 - crpix2) * cdelt2

        # convert brightness units to Jy/pixel
        if 'JY/PIXEL' in bunit.upper():
            pass
        else:
            print 'cannot handle BUNIT=%s' % bunit

        hdulist.close()

        # pixel area in sr
        fits_pixarea = abs(np.deg2rad(cdelt1 / 3600) * 
          np.deg2rad(cdelt2 / 3600)) 

        print 'magnifying source by factor', self.magnification
        magnified_axis1 = fits_axis1 * self.magnification
        magnified_axis2 = fits_axis2 * self.magnification

        # calculate pixel size of desired cube, should be as small as
        # the Nyquist freq of the longest baseline / highest freq
        target_pixsize = 1.0 / (self.wn_max * 100.0 * self.max_baseline)
        target_pixsize /= 2.0
        target_pixarea = target_pixsize**2
        target_pixsize = np.rad2deg(target_pixsize) * 3600
        print 'target pixel size', target_pixsize, 'arcsec'

        # now get size of target image to be created. This will be 
        # square, big enough to cover the primary beam, and centred at 
        # the FITS data array centre
        primary_beam_r = 1.22 / (self.wn_min * 100 * self.m1_diameter)
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

        smooth = abs(target_pixsize / (cdelt1 * self.magnification))
        print 'smoothing parameter', smooth

        # interpolate from FITS image onto target image
        # ..RectBivariateSpline requires axis1 and axis2 be monotonically
        # ..increasing, flip if need to
        axis1_mult = np.sign(magnified_axis1[1] - magnified_axis1[0])
        magnified_axis1 *= axis1_mult

        axis2_mult = np.sign(magnified_axis2[1] - magnified_axis2[0])
        magnified_axis2 *= axis2_mult

        print np.shape(magnified_axis1), np.shape(magnified_axis2), np.shape(fits_data)
        print '..setting up the interpolation'
        interp = interpolate.RectBivariateSpline(magnified_axis2, 
          magnified_axis1, fits_data, s=smooth)
        print '..interpolating'
        target_data = interp(target_axis, target_axis, grid=True)

        # multiply interpolated image by ratio of pixel areas to make
        # units brightness/pixel for the new pixel size
        print target_pixarea / (fits_pixarea * self.magnification**2)

        target_data *= (target_pixarea / (fits_pixarea * self.magnification**2))

        if self.plotname:
            # plot the image
            print '..plotting'
            plt.figure()

            plt.subplot(211)
            plt.imshow(fits_data, interpolation='nearest', origin='lower',
              aspect='equal', extent=[fits_axis1[0], fits_axis1[-1],
              fits_axis2[0], fits_axis2[-1]],
              vmax=np.max(fits_data)/5,
              vmin=np.min(fits_data))
            plt.colorbar(orientation='vertical')
            plt.axis('image')
            plt.title(object_name)

            plt.subplot(212)
            plt.imshow(target_data, interpolation='nearest', origin='lower',
              aspect='equal', extent=[target_axis[0], target_axis[-1],
              target_axis[0], target_axis[-1]],
              vmax=np.max(target_data)/5, vmin=np.min(target_data))
            plt.colorbar(orientation='vertical')
            plt.axis('image')
            plt.title('%s - interpolated' % object_name)

            filename = self.plotname
            plt.savefig(filename)
            plt.close()

        axis1 = co.Axis(data=target_axis, title=ctype1, units='arcsec') 
        axis2 = co.Axis(data=target_axis, title=ctype2, units='arcsec') 
        image = co.Image(target_data, title=object_name, axes=[axis1, axis2],
          units='Jy/pixel', info={'wn nominal':self.wn_nominal})
        return image


class MakeModelThinRing(object):
    """
    """

    def __init__(self, name, radius, tilt, xcentre, ycentre, max_baseline,
     wn_min, wn_max,  m1_diameter=2.0, verbose=False, plotname=None):
        self.name = name
        self.radius = radius
        self.tilt = tilt
        self.xcentre = xcentre
        self.ycentre = ycentre
        self.max_baseline = max_baseline
        self.wn_max = wn_max
        self.wn_min = wn_min
        self.m1_diameter = m1_diameter
        self.verbose = verbose
        self.plotname = plotname

    def run(self):
        print 'MakeModelImage.run'

        # calculate pixel size of desired cube, should be as small as
        # the Nyquist freq of the longest baseline / highest freq
        pixsize = 1.0 / (self.wn_max * 100.0 * self.max_baseline)
        pixsize /= 2.0
        pixsize = np.rad2deg(pixsize) * 3600
        print 'target pixel size', pixsize, 'arcsec'

        # now get size of target image to be created. This will be 
        # square, big enough to cover the primary beam, and centred at 
        # the FITS data array centre
        primary_beam_r = 1.22 / (self.wn_min * 100 * self.m1_diameter)
        primary_beam_r = np.rad2deg(primary_beam_r) * 3600 
        target_extent = 2 * primary_beam_r
        print 'primary beam diameter', primary_beam_r, 'arcsec'

        tnpix = int(target_extent / pixsize)

        trpix = int(tnpix / 2)
        trval = 0
        tdelt = pixsize
        target_axis = trval + (np.arange(tnpix, dtype=np.float) - trpix) * \
          tdelt 

        print 'target image will be', tnpix, 'x', tnpix, 'square'

        target_data = np.zeros([tnpix, tnpix])

        theta = np.arange(100) * 2 * np.pi / 100
        xlist = self.radius * np.sin(theta) 
        ylist = self.radius * np.cos(theta)
        ylist *= np.cos(np.deg2rad(self.tilt))

        xlist += self.xcentre
        ylist += self.ycentre
         
        for i,x in enumerate(xlist):
            y = ylist[i]

            ypix = np.argmin(np.abs(target_axis - y))
            xpix = np.argmin(np.abs(target_axis - x))

            ylo = int(np.floor(ypix - 4))
            yhi = int(np.ceil(ypix + 4))
            xlo = int(np.floor(xpix - 4))
            xhi = int(np.ceil(xpix + 4))

            for ypix in range(ylo, yhi):
                for xpix in range(xlo, xhi):
                    if ypix < 0 or ypix > len(target_axis)-1 or \
                      xpix < 0 or xpix > len(target_axis)-1:
                        continue
                    r2 = (target_axis[ypix] - y)**2 + (target_axis[xpix] - x)**2
                    r2 /= (target_axis[1] - target_axis[0])**2
                    target_data[ypix, xpix] += np.exp(-r2)
#              np.abs(target_axis - ylist[i]) < tdelt/2,
#              np.abs(target_axis - x) < tdelt/2] += 1.0

        if self.plotname is not None:

            # plot the image
            plt.figure()

            plt.imshow(target_data, interpolation='nearest', origin='lower',
              aspect='equal', extent=[target_axis[0], target_axis[-1],
              target_axis[0], target_axis[-1]],
              vmax=np.max(target_data), vmin=np.min(target_data))
            plt.colorbar(orientation='vertical')
            plt.axis('image')
            plt.title(self.name)

            plt.savefig(self.plotname)
            plt.close()

        axis1 = co.Axis(data=target_axis, title='RA', units='arcsec') 
        axis2 = co.Axis(data=target_axis, title='Dec', units='arcsec') 
        image = co.Image(target_data, title=self.name, axes=[axis2, axis1])
        return image


class MakeModelThickRing(object):
    """
    """

    def __init__(self, name, rinner, router, tilt, rot, xcentre, ycentre,
     max_baseline, wn_min, wn_max,  m1_diameter=2.0, verbose=False,
     plotname=None):
        self.name = name
        self.rinner = rinner
        self.router = router
        self.tilt = tilt
        self.rot = rot
        self.xcentre = xcentre
        self.ycentre = ycentre
        self.max_baseline = max_baseline
        self.wn_max = wn_max
        self.wn_min = wn_min
        self.m1_diameter = m1_diameter
        self.verbose = verbose
        self.plotname = plotname

    def run(self):
        print 'MakeModelThickRing.run'

        # calculate pixel size of desired cube, should be as small as
        # the Nyquist freq of the longest baseline / highest freq
        pixsize = 1.0 / (self.wn_max * 100.0 * self.max_baseline)
        pixsize /= 2.0
        pixsize = np.rad2deg(pixsize) * 3600
        print 'target pixel size', pixsize, 'arcsec'

        # now get size of target image to be created. This will be 
        # square, big enough to cover the primary beam, and centred at 
        # the FITS data array centre
        print self.wn_min, self.m1_diameter
        primary_beam_r = 1.22 / (self.wn_min * 100 * self.m1_diameter)
        primary_beam_r = np.rad2deg(primary_beam_r) * 3600 
        target_extent = 2 * primary_beam_r
        print 'primary beam diameter', primary_beam_r, 'arcsec'

        tnpix = int(target_extent / pixsize)

        trpix = int(tnpix / 2)
        trval = 0
        tdelt = pixsize
        target_axis = trval + (np.arange(tnpix, dtype=np.float) - trpix) * \
          tdelt 

        print 'target image will be', tnpix, 'x', tnpix, 'square'

        target_data = np.zeros([tnpix, tnpix])

        # now make the thick ring
        theta = np.arange(200) * 2 * np.pi / 200
        rstep = (self.router - self.rinner) / 100

        xlist = []
        ylist = []
        for ir in np.arange(100):
            r = self.rinner + rstep * ir

            xlist += list(r * np.sin(theta)) 
            ylist += list(r * np.cos(theta))

        # projection on sky
        xlist = np.array(xlist)
        ylist = np.array(ylist)

        # tilt
        ylist *= np.cos(np.deg2rad(self.tilt))

        # rotation
        yrlist = ylist * np.cos(np.deg2rad(self.rot)) - xlist * np.sin(np.deg2rad(self.rot))
        xrlist = ylist * np.sin(np.deg2rad(self.rot)) + xlist * np.cos(np.deg2rad(self.rot))

        # offset from centre
        xrlist += self.xcentre
        yrlist += self.ycentre
         
        for i,x in enumerate(xrlist):
            y = yrlist[i]

            ypix = np.argmin(np.abs(target_axis - y))
            xpix = np.argmin(np.abs(target_axis - x))

            ylo = int(np.floor(ypix - 4))
            yhi = int(np.ceil(ypix + 4))
            xlo = int(np.floor(xpix - 4))
            xhi = int(np.ceil(xpix + 4))

            for ypix in range(ylo, yhi):
                for xpix in range(xlo, xhi):
                    if ypix < 0 or ypix > len(target_axis)-1 or \
                      xpix < 0 or xpix > len(target_axis)-1:
                        continue
                    r2 = (target_axis[ypix] - y)**2 + (target_axis[xpix] - x)**2
                    r2 /= (target_axis[1] - target_axis[0])**2
                    target_data[ypix, xpix] += np.exp(-r2)

        if self.plotname is not None:

            # plot the image
            plt.figure()

            plt.imshow(target_data, interpolation='nearest', origin='lower',
              aspect='equal', extent=[target_axis[0], target_axis[-1],
              target_axis[0], target_axis[-1]],
              vmax=np.max(target_data), vmin=np.min(target_data))
            plt.colorbar(orientation='vertical')
            plt.axis('image')
            plt.title(self.name)

            plt.savefig(self.plotname)
            plt.close()

        axis1 = co.Axis(data=target_axis, title='RA', units='arcsec') 
        axis2 = co.Axis(data=target_axis, title='Dec', units='arcsec') 
        image = co.Image(target_data, title=self.name, axes=[axis2, axis1])
        return image


class MakeModelComet(object):
    """
    """

    def __init__(self, name, xcentre, ycentre,
     max_baseline, wn_min, wn_max,  m1_diameter=2.0, verbose=False,
     plotname=None):
        self.name = name
        self.xcentre = xcentre
        self.ycentre = ycentre
        self.max_baseline = max_baseline
        self.wn_max = wn_max
        self.wn_min = wn_min
        self.m1_diameter = m1_diameter
        self.verbose = verbose
        self.plotname = plotname

    def run(self):
        print 'MakeModelComet.run'

        # calculate pixel size of desired cube, should be as small as
        # the Nyquist freq of the longest baseline / highest freq
        pixsize = 1.0 / (self.wn_max * 100.0 * self.max_baseline)
        pixsize /= 2.0
        pixsize = np.rad2deg(pixsize) * 3600
        print 'target pixel size', pixsize, 'arcsec'

        # now get size of target image to be created. This will be 
        # square, big enough to cover the primary beam, and centred at 
        # the FITS data array centre
        print self.wn_min, self.m1_diameter
        primary_beam_r = 1.22 / (self.wn_min * 100 * self.m1_diameter)
        primary_beam_r = np.rad2deg(primary_beam_r) * 3600 
        target_extent = 2 * primary_beam_r
        print 'primary beam diameter', primary_beam_r, 'arcsec'

        tnpix = int(target_extent / pixsize)

        trpix = int(tnpix / 2)
        trval = 0
        tdelt = pixsize
        target_axis = trval + (np.arange(tnpix, dtype=np.float) - trpix) * \
          tdelt 

        print 'target image will be', tnpix, 'x', tnpix, 'square'

        target_data = np.zeros([tnpix, tnpix])

        # now make the comet
        thetas = (np.arange(200) * np.pi / 200) - (np.pi / 2)

        i_centre = np.argmin(np.abs(target_axis - self.xcentre))
        j_centre = np.argmin(np.abs(target_axis - self.ycentre))
        k_centre = tnpix / 2
        cube = np.zeros([tnpix,tnpix,tnpix])

        jets = [(np.pi/3, 0.0, 100), (-4*np.pi/5, np.pi/5, 50)]

        for k in np.arange(tnpix, dtype=np.float):
            print k
            for j in np.arange(tnpix, dtype=np.float):
                i = np.arange(tnpix, dtype=np.float)
                r2 = (i-i_centre)**2 + (j-j_centre)**2 + (k-k_centre)**2
                rij = np.sqrt((i-i_centre)**2 + (j-j_centre)**2)

                theta = np.arctan2((k-k_centre), rij)
                phi = np.arctan2((j-j_centre), (i-i_centre))

                density = np.ones(i.shape)
                density[r2>1] = 1.0 / r2[r2>1]

                for jet in jets:
                    phi_jet = jet[0]
                    theta_jet = jet[1]
                    d_jet = jet[2]

                    cos_dist_jet = np.sin(theta) * np.sin(theta_jet) + \
                      np.cos(theta) * np.cos(theta_jet) * np.cos(phi - phi_jet)
                    ind = np.logical_and(cos_dist_jet > np.cos(3 * np.pi / 180), (r2>1))
                    density[ind] += d_jet / r2[ind]

                cube[k,j,:] = density 

        target_data = np.sum(cube, axis=0) 

        if self.plotname is not None:

            # plot the image
            plt.figure()

            plt.imshow(target_data, interpolation='nearest', origin='lower',
              aspect='equal', extent=[target_axis[0], target_axis[-1],
              target_axis[0], target_axis[-1]],
              vmax=np.max(target_data), vmin=np.min(target_data))
            plt.colorbar(orientation='vertical')
            plt.axis('image')
            plt.title(self.name)

            plt.savefig(self.plotname)
            plt.close()

        axis1 = co.Axis(data=target_axis, title='RA', units='arcsec') 
        axis2 = co.Axis(data=target_axis, title='Dec', units='arcsec') 
        image = co.Image(target_data, title=self.name, axes=[axis2, axis1])
        return image


class MakeSpectrum(object):
    """
    """

    def __init__(self, temperature, beta, wn_min, wn_max, wn_step,
      peak_flux, spectral_features=[], plotname=None):
        assert wn_max > wn_min
        assert wn_step > 0

        self.temperature = temperature
        self.beta = beta
        self.wn_min = wn_min
        self.wn_max = wn_max
        self.wn_step = wn_step
        self.peak_flux = peak_flux
        self.plotname = plotname

        self.forsterite = [
          (67,0.11),
          (67.1,0.075),
          (67.2,0.06),
          (67.3,0.13),
          (67.4,0),
          (67.5,-0.05),
          (67.6,0.17),
          (67.7,0.35),
          (67.8,0.275),
          (67.9,0.2),
          (68.0,0.15),
          (68.1,0.22),
          (68.2,0.125),
          (68.3,0.05),
          (68.4,0.1),
          (68.5,0.25),
          (68.6,0.225),
          (68.7,0.25),
          (68.8,0.5),
          (68.9,0.9),
          (69.0,1.2),
          (69.1,1.38),
          (69.2,1.5),
          (69.3,1.55),
          (69.4,1.5),
          (69.5,1.22),
          (69.6,1.15),
          (69.7,0.95),
          (69.8,0.7),
          (69.9,0.6),
          (70.0,0.55),
          (70.1,0.45),
          (70.2,0.4),
          (70.3,0.45),
          (70.4,0.6),
          (70.5,0.65),
          (70.6,0.44),
          (70.7,0.2),
          (70.8,0.04),
          (70.9,-0.1),
          (71.0,-0.15)]

        self.protostar = [
          (63.0,0.1),
          (63.1,6.5),
          (63.2,15.0),
          (63.3,12.1),
          (63.4,0.1),
          (66.2,0.1),
          (66.3,1),
          (66.4,9.0),
          (66.5,1.5),
          (66.6,0.1),
          (71.8,0.1),
          (71.9,8.5),
          (72.0,2.0),
          (72.1,3.5),
          (72.2,2.0),
          (72.3,0.1),
          (100.3,0.1),
          (100.4,5),
          (100.5,40),
          (100.6,5),
          (100.7,0.1),
          (153.2,0.1),
          (153.3,10),
          (153.4,85.0),
          (153.5,10),
          (153.6,0.1),
          (324.8,0.1),
          (324.9,8),
          (325.0,80),
          (325.1,8),
          (325.2,0.1)]

        self.comet = [
#          (66.0, 5e-3),
#          (300.0, 5e-3)]
          (66.0, 400*5e-3),
          (300.0, 400*5e-3)]

        self.spectral_features = []
        for feature in spectral_features:
            if isinstance(feature, basestring):
                if 'FORSTERITE' in feature.upper():
                    self.spectral_features.append(self.forsterite)
                elif 'PROTOSTAR' in feature.upper():
                    self.spectral_features.append(self.protostar)
                elif 'COMET' in feature.upper():
                    self.spectral_features.append(self.comet)
                else:
                    raise Exception, 'unknown feature named'
            else:
                self.spectral_features.append(feature)

    def run(self):
        print 'MakeSpectrum.run'

        nwn = np.ceil(abs((self.wn_max - self.wn_min) / self.wn_step))
        wn_vector = self.wn_min + np.arange(nwn) * self.wn_step

        # spectrum
        bb = BB_spectrum(self.temperature, wn_vector, cutoffmin=wn_vector[2],
          cutoffmax=wn_vector[-2], emissivity=1.0)
        spectrum = bb.calculate()

        # now multiply by (lambda / lambda_normalise)**beta to mimic
        # dust emissivity variation
        emissivity = np.ones(np.shape(wn_vector))
        emissivity *= (wn_vector[0] / wn_vector) ** self.beta
        spectrum *= emissivity

        # normalise the spectrum to the given peak
        norm = self.peak_flux / np.max(spectrum)
        spectrum *= norm

        # now add in spectral features
        for feature in self.spectral_features:
            lambs = []
            fluxes = []
            for (lamb,flux) in feature:
                lambs.append(lamb)
                fluxes.append(flux)
            lambs = np.array(lambs)
            fluxes = np.array(fluxes)

            # assume lambda is in microns
            wn_width = np.abs(wn_vector[1] - wn_vector[0])
            for iwn,wn in enumerate(wn_vector):
                wn_lo = wn - wn_width / 2
                wn_hi = wn + wn_width / 2
                lam_lo = 1e4 / wn_hi
                lam_hi = 1e4 / wn_lo

                feat = fluxes[np.logical_and(lambs > lam_lo, lambs < lam_hi)]
                print feat
                if len(feat):
                    feat = np.sum(feat) / len(feat)
                    spectrum[iwn] += feat

        # bb value at wn_normalise
#        bb = BB_spectrum(self.temperature, [self.wn_normalise], emissivity=1.0)
#        norm = bb.calculate()
#        norm = self.flux_normalise / norm

        # and normalise the spectrum at wn_normalise
#        spectrum *= norm

        if self.plotname is not None:
            # plot the spectrum
            plt.figure()

            plt.plot(wn_vector, spectrum, 'k')
            plt.savefig(self.plotname)
            plt.close()

        axis = co.Axis(data=wn_vector, title='wavenumber', units='cm-1')
        result = co.Spectrum(spectrum, axis=axis, title='')
        return result
