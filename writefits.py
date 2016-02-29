"""This module contains classes and methods used to save data from the
simulated observation to a FITS file. 
"""

from __future__ import absolute_import

import numpy as np
import os.path
import astropy.io.fits as pyfits


class WriteFITS(object):
    """Class to write out observed timeline in a FITS Table.

    Contains methods:
    __init__
    run
    __repr__
    """

    def __init__(self, previous_results):
        """Constructor.

        Keyword parameters:
        previous_results - simulator result structure.
        """
        self.previous_results = previous_results
        self.result = {}

    def run(self):
        """Method to write the FITS data.
        """
        #print 'WriteFITS.run'

        # construct the name of the file
        runtime = self.previous_results['runtime']
        fitsname = '%s.fits' % runtime

        # get specific items from the results that will be need in
        # the reduction
        fts = self.previous_results['fts']
        smec_opd_to_mpd = fts['smec_opd_to_mpd']
        wnmin = fts['wnmin']

        skymodel = self.previous_results['skymodel']
        spatial_axis = skymodel['spatial axis [arcsec]']
        npix = len(spatial_axis)

        # get list of instrument observations
        observe = self.previous_results['observe']
        obs_timeline = observe['observed_timeline']
        observed_times = obs_timeline.keys()
        observed_times.sort()

        # construct lists of the values to be stored in each Table column
        timelist = []
        scan = []
        baseline_x = []
        baseline_y = []
        baseline_z = []
        smec_position = []
        smec_nominal_position = []
        smec_velocity_error = []
        flag = []
        data = []
        pure_data = []
        cr_data = []
        detector_noise_data = []
        pointing1_x = []
        pointing1_y = []
        pointing2_x = []
        pointing2_y = []

        for t in observed_times:
            config = obs_timeline[t]

            timelist.append(config.time)
            scan.append(config.scan_number)
            baseline_x.append(config.baseline_x)
            baseline_y.append(config.baseline_y)
            baseline_z.append(config.baseline_z)
            smec_position.append(config.smec_position)
            smec_nominal_position.append(config.smec_nominal_position)
            smec_velocity_error.append(config.smec_vel_error)
            flag.append(config.flag)
            data.append(config.data)
            pure_data.append(config.pure_data)
            cr_data.append(config.cr_data)
            detector_noise_data.append(config.detector_noise_data)
            pointing1_x.append(config.pointing1_x)
            pointing1_y.append(config.pointing1_y)
            pointing2_x.append(config.pointing2_x)
            pointing2_y.append(config.pointing2_y)

        # create a Header object and primary HDU - this will contain
        # general information
        prihdr = pyfits.Header()
        prihdr['COMMENT'] = 'This FITS file was created by pyfiins at %s' % \
          runtime
        prihdr['DATE'] = runtime
        prihdr['smec_o2m'] = smec_opd_to_mpd 
        prihdr['wnmin'] = wnmin
        prihdr['npix'] = npix

        prihdu = pyfits.PrimaryHDU(header=prihdr)

        # create list of Header Data Unit objects, include the primary HDU
        hdulist = pyfits.HDUList([prihdu])

        # create an ImageHDU to contain the spatail axis array - might be a
        # better way to do this but I can't find it
        spatial_hdu = pyfits.Header()
        spatial_hdu['comment'] = \
          'This data unit contains the spatial axis information'
        spatial_hdu = pyfits.ImageHDU(data=spatial_axis, header=spatial_hdu)
        hdulist.append(spatial_hdu)

        # create an HDU to contain the Table and append it to the list
        hdulist.append(pyfits.BinTableHDU.from_columns(
          pyfits.ColDefs([
          pyfits.Column(name='Time', format='D',
          array=np.array(timelist)),
          pyfits.Column(name='Scan', format='I',
          array=np.array(scan)),
          pyfits.Column(name='Baseline x', format='E',
          array=np.array(baseline_x)),
          pyfits.Column(name='Baseline y', format='E',
          array=np.array(baseline_y)),
          pyfits.Column(name='Baseline z', format='E',
          array=np.array(baseline_z)),
          pyfits.Column(name='SMEC Position', format='E',
          array=np.array(smec_position)),
          pyfits.Column(name='SMEC Nominal Position', format='E',
          array=np.array(smec_nominal_position)),
          pyfits.Column(name='SMEC Velocity Error', format='E',
          array=np.array(smec_velocity_error)),
          pyfits.Column(name='Flag', format='L',
          array=np.array(flag)),
          pyfits.Column(name='Data', format='E',
          array=np.array(data)),
          pyfits.Column(name='Pure Data', format='E',
          array=np.array(pure_data)),
          pyfits.Column(name='Cosmic Ray Data', format='E',
          array=np.array(cr_data)),
          pyfits.Column(name='Detector Noise Data', format='E',
          array=np.array(detector_noise_data)),
          pyfits.Column(name='Pointing1 X', format='E',
          array=np.array(pointing1_x)),
          pyfits.Column(name='Pointing1 Y', format='E',
          array=np.array(pointing1_y)),
          pyfits.Column(name='Pointing2 X', format='E',
          array=np.array(pointing2_x)),
          pyfits.Column(name='Pointing2 Y', format='E',
          array=np.array(pointing2_y))])))

        # write the HDU list to a file
        hdulist.writeto(fitsname, clobber=True)
        self.result['fitsfile'] = fitsname

        return self.result

    def __repr__(self):
        return '''
WriteFITS:
  FITS file : {fitsfile}
'''.format(fitsfile=self.result['fitsfile'])


class WriteFITSCube(object):
    """Class to write cube to a FITS primary HDU.

    Contains methods:
    __init__
    run 
    """

    def __init__(self, cube, fitsname='writefitscube.fits'):
        """Constructor.

        Keyword parameters:
        cube     - the cube to be stored
        fitsname - name of FITS file (default 'writefitscube.fits')
        """
        self.cube = cube
        self.fitsname = fitsname
        self.result = {}

    def run(self):
        """Method to write the FITS data.
        """
        print 'WriteFITSCube.run'

        # create a Header object and primary HDU - this will contain
        # general information
        prihdr = pyfits.Header()
        prihdr['COMMENT'] = 'This FITS file was created by WriteFITSCube'

        # data cube comes in with indices [wn, dec, ra]. FITS file expects
        # data in Fortran order [ra, dec, wn].
        # Do some axis swapping to achieve this.
        data = self.cube.data

        axis3 = self.cube.axes[0].data
        title3 = self.cube.axes[0].title
        units3 = self.cube.axes[0].units
        axis2 = self.cube.axes[1].data
        title2 = self.cube.axes[1].title
        units2 = self.cube.axes[1].units
        axis1 = self.cube.axes[2].data
        title1 = self.cube.axes[2].title
        units1 = self.cube.axes[2].units

        # calculate some header values
        crpix1 = float(1)
        crpix2 = float(1)
        crpix3 = float(1)

        if 'ARCSEC' in units1.upper():
            cdelt1 = (axis1[1] - axis1[0]) / 3600.0
            crval1 = axis1[crpix1-1] / 3600.0
            cunit1 = 'DEG     '
            ctype1 = 'RA---SIN'
        elif 'CM-1' in  units1.upper():
            cdelt1 = (axis1[1] - axis1[0]) * 3.0e10
            crval1 = axis1[crpix1-1] * 3.0e10
            cunit1 = 'HZ      '
            ctype1 = 'FREQ    '
        elif 'HZ' in  units1.upper():
            cdelt1 = (axis1[1] - axis1[0])
            crval1 = axis1[crpix1-1]
            cunit1 = 'HZ      '
            ctype1 = 'FREQ    '
        elif 'AU' in units1.upper():
            cdelt1 = axis1[1] - axis1[0]
            crval1 = axis1[crpix1-1]
            cunit1 = 'AU      '
            ctype1 = 'X-OFFSET'
        else:
            cdelt1 = axis1[1] - axis1[0]
            crval1 = axis1[crpix1-1]
            cunit1 = 'UNKNOWN '
            ctype1 = 'UNKNOWN '

        if 'ARCSEC' in units2.upper():
            cdelt2 = (axis2[1] - axis2[0]) / 3600.0
            crval2 = axis2[crpix2-1] / 3600.0
            cunit2 = 'DEG     '
            ctype2 = 'DEC--SIN'
        elif 'CM-1' in  units2.upper():
            cdelt2 = (axis2[1] - axis2[0]) * 3.0e10
            crval2 = axis2[crpix2-1] * 3.0e10
            cunit2 = 'HZ      '
            ctype2 = 'FREQ    '
        elif 'HZ' in  units1.upper():
            cdelt2 = (axis2[1] - axis2[0])
            crval2 = axis2[crpix2-1]
            cunit2 = 'HZ      '
            ctype2 = 'FREQ    '
        elif 'AU' in units2.upper():
            cdelt2 = axis2[1] - axis2[0]
            crval2 = axis2[crpix2-1]
            cunit2 = 'AU      '
            ctype2 = 'Y-OFFSET'
        else:
            cdelt2 = axis2[1] - axis2[0]
            crval2 = axis2[crpix2-1]
            cunit2 = 'UNKNOWN '
            ctype2 = 'UNKNOWN '

        if 'ARCSEC' in units3.upper():
            cdelt3 = (axis3[1] - axis3[0]) / 3600.0
            crval3 = axis3[crpix3-1] / 3600.0
            cunit3 = 'DEG     '
            ctype3 = 'UNKNOWN '
        elif 'CM-1' in  units3.upper():
            cdelt3 = (axis3[1] - axis3[0]) * 3.0e10
            crval3 = axis3[crpix3-1] * 3.0e10
            cunit3 = 'HZ      '
            ctype3 = 'FREQ    '
        elif 'HZ' in  units3.upper():
            cdelt3 = (axis3[1] - axis3[0])
            crval3 = axis3[crpix3-1]
            cunit3 = 'HZ      '
            ctype3 = 'FREQ    '
        elif 'AU' in units3.upper():
            cdelt3 = axis3[1] - axis3[0]
            crval3 = axis3[crpix3-1]
            cunit3 = 'AU      '
            ctype3 = 'OFFSET  '
        else:
            cdelt3 = axis3[1] - axis3[0]
            crval3 = axis3[crpix3-1]
            cunit3 = 'UNKNOWN '
            ctype3 = 'UNKNOWN '

        prihdr['CDELT1'] = cdelt1
        prihdr['CDELT2'] = cdelt2
        prihdr['CDELT3'] = cdelt3

        prihdr['CRVAL1'] = crval1
        prihdr['CRVAL2'] = crval2
        prihdr['CRVAL3'] = crval3

        prihdr['CRPIX1'] = crpix1
        prihdr['CRPIX2'] = crpix2
        prihdr['CRPIX3'] = crpix3

        prihdr['CUNIT1'] = cunit1
        prihdr['CUNIT2'] = cunit2
        prihdr['CUNIT3'] = cunit3

        prihdr['CTYPE1'] = ctype1
        prihdr['CTYPE2'] = ctype2
        prihdr['CTYPE3'] = ctype3
   
        prihdr['BTYPE'] = 'Intensity'
        prihdr['BUNIT'] = 'JY/PIXEL'

        # some fake coordinate info
        prihdr['EQUINOX'] = float(2000)
        prihdr['RADESYS'] = 'FK5     '
        prihdr['SPECSYS'] = 'LSRK    '
        # VELREF==1 implies LSR velocity frame in CASA
        prihdr['VELREF'] = 1

        prihdu = pyfits.PrimaryHDU(header=prihdr, data=data)
        hdulist = pyfits.HDUList([prihdu])

        # write the HDU list to a file
        hdulist.writeto(self.fitsname, clobber=True)
        self.result['fitsfile'] = self.fitsname

        return self.result
