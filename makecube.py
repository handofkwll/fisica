from __future__ import absolute_import

import numpy as np
import scipy.interpolate as interpolate
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt

import common.commonobjects as co
from skygenerator import BB_spectrum
import writefits


class MakeCube(object):
    """Class to make a cube from a 2-d image crossed with a spectrum. 
    """

    def __init__(self, image, spectrum):
        self.image = image
        self.spectrum = spectrum

    def run(self):
        print 'MakeCube.run'

        # the cube shape is just the image replicated through the spectral
        # dimension
        cubeshape = self.image.data.shape + self.spectrum.data.shape
        cube = np.zeros(cubeshape, np.float)
      
        # populate each plane of the cube in turn
        for iwn,wn in enumerate(self.spectrum.axis.data):
            cube[:,:,iwn] = self.image.data * self.spectrum.data[iwn]

        # write the cube to FITS
        cube = co.Image(data=cube, axes=self.image.axes + [self.spectrum.axis])
        cubewriter = writefits.WriteFITSCube(cube, 'cube.fits')
        cubewriter.run() 
        

class MakeImage(object):
    """
    """

    def __init__(self, fits_image, max_baseline, wn_nominal, wn_max, wn_min,
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


class MakeModelImage(object):
    """
    """

    def __init__(self, name, radius, tilt, xcentre, ycentre, max_baseline,
     wn_max, wn_min,  m1_diameter=2.0, verbose=False):
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
            target_data[np.abs(target_axis - x) < tdelt/2,
                   np.abs(target_axis - ylist[i]) < tdelt/2] += 1.0

        # plot the image
        plt.figure()

        plt.imshow(target_data, interpolation='nearest', origin='lower',
          aspect='equal', extent=[target_axis[0], target_axis[-1],
          target_axis[0], target_axis[-1]],
          vmax=np.max(target_data), vmin=np.min(target_data))
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title(self.name)

        filename = 'makeimage.png'
        plt.savefig(filename)
        plt.close()

        axis1 = co.Axis(data=target_axis, title='RA', units='arcsec') 
        axis2 = co.Axis(data=target_axis, title='Dec', units='arcsec') 
        image = co.Image(target_data, title=self.name, axes=[axis1, axis2])
        return image


class MakeSpectrum(object):
    """
    """

    def __init__(self, temperature, beta, wn_min, wn_max, wn_step,
      wn_normalise, plotname=None):
        assert wn_max > wn_min
        assert wn_step > 0

        self.temperature = temperature
        self.beta = beta
        self.wn_min = wn_min
        self.wn_max = wn_max
        self.wn_step = wn_step
        self.wn_normalise = wn_normalise
        self.plotname = plotname

    def run(self):
        print 'MakeSpectrum.run'

        nwn = np.ceil(abs((self.wn_max - self.wn_min) / self.wn_step))
        wn_vector = self.wn_min + np.arange(nwn) * self.wn_step

        # spectrum
        bb = BB_spectrum(self.temperature, wn_vector, cutoffmin=wn_vector[2],
          cutoffmax=wn_vector[-2], emissivity=1.0)
        spectrum = bb.calculate()

        # bb value at wn_normalise
        bb = BB_spectrum(self.temperature, [self.wn_normalise], emissivity=1.0)
        norm = bb.calculate()

        # and normalise the spectrum at wn_normalise
        spectrum /= norm

        # now multiply by (lambda / lambda_normalise)**beta to mimic
        # dust emissivity variation
        emissivity = np.ones(np.shape(wn_vector))
        emissivity *= (self.wn_normalise / wn_vector) ** self.beta
        spectrum *= emissivity

        if self.plotname is not None:
            # plot the spectrum
            plt.figure()

            plt.plot(wn_vector, spectrum, 'k')
            plt.plot(wn_vector, spectrum / emissivity, 'b')
            plt.savefig(self.plotname)
            plt.close()

        axis = co.Axis(data=wn_vector, title='wavenumber', units='cm-1')
        result = co.Spectrum(spectrum, axis=axis, title='')
        return result
