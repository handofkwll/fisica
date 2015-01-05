from __future__ import absolute_import

import collections
import math
import numpy as np
import pp

import common.commonobjects as co

def dirty(b_x_list, b_y_list, spectra, wn_spectra, iwn, wn, spatial_axis, npix):
    image = numpy.zeros([npix, npix], numpy.float)
    beam = numpy.zeros([npix, npix], numpy.float)

    for ibx,b_x in enumerate(b_x_list):
        b_x = b_x_list[ibx]
        b_y = b_y_list[ibx]
        spectrum = spectra[ibx]
        wn_spectrum = wn_spectra[ibx]

        # must be a better way of doing this
        argx = numpy.radians(b_x * (wn * 100.0) * spatial_axis / 3600.0)
        argy = numpy.radians(b_y * (wn * 100.0) * spatial_axis / 3600.0)
        valx = numpy.exp(2.0j * math.pi * argx)
        valy = numpy.exp(2.0j * math.pi * argy)

        vis = spectrum[wn_spectrum==wn]

        # numpy arrays are [row,col] or here [y,x]
        for iy,y in enumerate(spatial_axis):
            contribution = (vis * valx * valy[iy]) + \
              numpy.conjugate(vis * valx * valy[iy])
            image[iy,:] += contribution.real
            contribution = (valx * valy[iy]) + \
              numpy.conjugate(valx * valy[iy])
            beam[iy,:] += contribution.real
   
    # normalise
    # The sum of the dirty beam should equal the 0 baseline vis measurement: 0
    # Normalising the volume would preserve flux under convolution but
    # is problematic as the volume is 0, so normalise the peak to 1.
    beam /= 2.0 * len(b_x_list)

    # likewise for dirty map
    image /= 2.0 * len(b_x_list)

    return image, beam


class DirtyImage(object):
    """Class to compute dirty image cube from uvspectra.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        print 'DirtyImage.run'

        # get relevant fts info
        fts = self.previous_results['fts']
        fts_wnmin = fts['wnmin']

        # get primary beam info
        cubeparameters = self.previous_results['cubeparameters']
        npix = cubeparameters['npix']
        spatial_axis = cubeparameters['spatial axis [arcsec]']

        # get observation list
        uvspectra = self.previous_results['uvspectra']['uvspectra']
        wavenumber = uvspectra[0].wavenumber
        wavenumber = wavenumber[wavenumber>fts_wnmin]

        # dirty image cube
        dirtyimage = np.zeros([npix, npix, len(wavenumber)], np.float)
        dirtybeam = np.zeros([npix, npix, len(wavenumber)], np.float)

        # calculate dirty image for each wn
        # uvspectrum objects don't pickle so unpack them first
        b_x_list = []
        b_y_list = []
        spectra = []
        wn_spectra = []
        for uvspectrum in uvspectra:
            b_x_list.append(uvspectrum.baseline_x)
            b_y_list.append(uvspectrum.baseline_y)
            spectra.append(uvspectrum.spectrum)
            wn_spectra.append(uvspectrum.wavenumber)

        jobs = {}
        for iwn,wn in enumerate(wavenumber):
            # submit jobs
            indata = (b_x_list, b_y_list, spectra, wn_spectra, iwn, wn, 
              spatial_axis, npix,)
            jobs[wn] = self.job_server.submit(dirty, indata, (),
              ('numpy','math',))

        for iwn,wn in enumerate(wavenumber):
            # collect and store results
            dirtyimage[:,:,iwn], dirtybeam[:,:,iwn] = jobs[wn]()

        self.result['dirtyimage'] = dirtyimage
        self.result['dirtybeam'] = dirtybeam
        self.result['spatial axis [arcsec]'] = spatial_axis
        self.result['spatial axis'] = spatial_axis
        self.result['wavenumber [cm-1]'] = wavenumber

        return self.result
                
    def __repr__(self):
        return '''
DirtyImage:
'''.format(
          num_uvspectra=len(self.uvspectra))

