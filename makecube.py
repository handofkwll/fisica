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

    def __init__(self, fits_image, max_baseline, wn_max, wn_min, magnification,
      m1_diameter=2.0, verbose=False):
        self.fits_image = fits_image
        self.max_baseline = max_baseline
        self.wn_max = wn_max
        self.wn_min = wn_min
        self.magnification = magnification
        self.m1_diameter = m1_diameter
        self.verbose = verbose

    def run(self):
        print 'MakeImage.run'

        hdulist = pyfits.open(self.fits_image)
        print hdulist.info()
        prihdr = hdulist[0].header
        if self.verbose:
            print repr(prihdr)
        image_data = hdulist[0].data

        # Note that, like C (and unlike FORTRAN), Python is 0-indexed and 
        # the indices have the slowest axis first and fastest changing axis
        # last; i.e. for a 2-D image, the fast axis (X-axis) which 
        # corresponds to the FITS NAXIS1 keyword, is the second index.

        nfreq = image_data.shape[0]

        # get rid of degenerate axes (freq), verify data are 2d and swap
        # axes - image_data ends with shape [AXIS2, AXIS1]
        image_data = np.squeeze(image_data)
        assert len(image_data.shape) == 2
        # don't swap axes - looks like Fomalhaut data is in C order
        #image_data = np.swapaxes(image_data, 0, 1)

        datashape = np.shape(image_data)
        nx = datashape[1]
        ny = datashape[0]
#        assert nx==ny

        try:
            object_name = prihdr['OBJECT']
        except:
            object_name = 'unknown'

        # pixel increments in degrees
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

        cunit3 = prihdr['CUNIT3']

        bunit = prihdr['BUNIT']

        # calculate axes of FITS image in arcsec
        if 'DEG' in cunit1.upper():
            # + 1 is in there because FITS arrays are 1-based
            axis1 = (np.arange(nx) + 1 - crpix1) * cdelt1 * 3600.0
        else:
            raise Exception, 'cannot handle CUNIT1=%s' % cunit1

        if 'DEG' in cunit2.upper():
            axis2 = (np.arange(ny) + 1 - crpix2) * cdelt2 * 3600.0
        else:
            raise Exception, 'cannot handle CUNIT2=%s' % cunit2

        # frequency axis
        if 'HZ' in cunit3.upper():
            # ..in Hz (+1 is in there because FITS arrays are 1-based)
            fvec = (np.arange(nfreq) + 1 - crpix3) * cdelt3 + crval3
            # ..in cm-1
            fvec_wn = fvec / 3.0e10
            # ..wavelength in metres
            fvec_wl = 3.0e8 / fvec
        elif 'M' in cunit3.upper():
            # ..wavelength in metres
            fvec_wl = (np.arange(nfreq) + 1 - crpix3) * cdelt3 + crval3
            # ..in Hz
            fvec = 3.0e8 / fvec_wl
            # ..in cm-1
            fvec_wn = fvec / 3.0e10
        else:
            raise Exception, 'cannot handle CUNIT3=%s' % cunit3

        # convert from Jy/pixel to W/m2/cm-1/sr-1
        if 'JY/PIXEL' in bunit.upper():
            pixsize_rad = np.deg2rad(abs(spatial_axis[1] - spatial_axis[0]) / 3600.0) 
        elif 'PW' in bunit.upper():
            print 'picoWatts'
        else:
            print 'cannot handle BUNIT=%s' % bunit

        hdulist.close()

        print 'magnifying source by factor', self.magnification
        axis1 *= self.magnification
        axis2 *= self.magnification

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

        # interpolate from FITS image onto target image
        # ..RectBivariateSpline requires axis1 and axis2 be monotonically
        # ..increasing, flip if need to
        axis1_mult = np.sign(axis1[1] - axis1[0])
        axis1 *= axis1_mult

        axis2_mult = np.sign(axis2[1] - axis2[0])
        axis2 *= axis2_mult

        interp = interpolate.RectBivariateSpline(axis1, axis2, image_data)
        target_data = interp(target_axis, target_axis, grid=True)

        # flip axes back again for plot
        axis1 *= axis1_mult
        axis2 *= axis2_mult

        # plot the image
        plt.figure()

        plt.subplot(211)
        plt.imshow(image_data, interpolation='nearest', origin='lower',
          aspect='equal', extent=[axis1[0], axis1[-1],
          axis2[0], axis2[-1]],
          vmax=np.max(image_data), vmin=np.min(image_data))
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title(object_name)

        plt.subplot(212)
        plt.imshow(target_data, interpolation='nearest', origin='lower',
          aspect='equal', extent=[target_axis[0], target_axis[-1],
          target_axis[0], target_axis[-1]],
          vmax=np.max(target_data), vmin=np.min(target_data))
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('%s - interpolated' % object_name)

        filename = 'makeimage.png'
        plt.savefig(filename)
        plt.close()

        axis1 = co.Axis(data=target_axis, title=ctype1, units='arcsec') 
        axis2 = co.Axis(data=target_axis, title=ctype2, units='arcsec') 
        image = co.Image(target_data, title=object_name, axes=[axis1, axis2])
        return image


class MakeSpectrum(object):
    """
    """

    def __init__(self, temperature, beta, wn_min, wn_max, wn_step):
        assert wn_max > wn_min
        assert wn_step > 0

        self.temperature = temperature
        self.beta = beta
        self.wn_min = wn_min
        self.wn_max = wn_max
        self.wn_step = wn_step

    def run(self):
        print 'MakeSpectrum.run'

        nwn = np.ceil(abs((self.wn_max - self.wn_min) / self.wn_step))
        wn_list = self.wn_min + np.arange(nwn) * self.wn_step

        print wn_list

        bb = BB_spectrum(self.temperature, wn_list, wn_list[2],
          wn_list[-2], 1.0)
        spectrum = bb.calculate()

        print spectrum

        # plot the image
        plt.figure()

        plt.plot(wn_list, spectrum)

        filename = 'makespectrum.png'
        plt.savefig(filename)
        plt.close()

        axis = co.Axis(data=wn_list, title='wavenumber', units='cm-1')
        result = co.Spectrum(spectrum, axis=axis, title='')
        return result
