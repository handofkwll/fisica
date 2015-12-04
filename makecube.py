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
         
        indices = np.arange(len(target_axis))
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


class MakeSpectrum(object):
    """
    """

    def __init__(self, temperature, beta, wn_min, wn_max, wn_step,
      flux_normalise, wn_normalise, spectral_features=[], plotname=None):
        assert wn_max > wn_min
        assert wn_step > 0

        self.temperature = temperature
        self.beta = beta
        self.wn_min = wn_min
        self.wn_max = wn_max
        self.wn_step = wn_step
        self.flux_normalise = flux_normalise
        self.wn_normalise = wn_normalise
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

        self.spectral_features = []
        for feature in spectral_features:
            if isinstance(feature, basestring) and \
              ('FORSTERITE' in feature.upper()):
                self.spectral_features.append(self.forsterite)
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

        # bb value at wn_normalise
        bb = BB_spectrum(self.temperature, [self.wn_normalise], emissivity=1.0)
        norm = bb.calculate()
        norm = self.flux_normalise / norm

        # and normalise the spectrum at wn_normalise
        spectrum *= norm

        # now multiply by (lambda / lambda_normalise)**beta to mimic
        # dust emissivity variation
        emissivity = np.ones(np.shape(wn_vector))
        emissivity *= (self.wn_normalise / wn_vector) ** self.beta
        spectrum *= emissivity

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
                if len(feat):
                    feat = np.sum(feat) / len(feat)
                    spectrum[iwn] += feat

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
