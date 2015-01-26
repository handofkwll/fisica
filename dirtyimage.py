from __future__ import absolute_import

import collections
import math
import numpy as np

def calculate_dirty_plane(b_x_list, b_y_list, spectra, wn_spectra,
  iwn, wn, spatial_axis, npix):
    """Routine to calculate one plane of a dirty image cube.

    Parameters:
    b_x_list     -
    b_y_list     -
    spectra      -
    wn_spectra   -
    iwn          - 
    wn           -
    spatial_axis -
    npix         - image plane has dimnesions [npix, npix]

    Returns:
    image        - [npix, npix] float array with dirty image 
    beam         - [npix, npix] float array with dirty beam

    """
    # arrays to hold the results
    image = numpy.zeros([npix, npix], numpy.float)
    beam = numpy.zeros([npix, npix], numpy.float)

    # iterate through the measured baselines
    for ibx,b_x in enumerate(b_x_list):
        b_x = b_x_list[ibx]
        b_y = b_y_list[ibx]
        spectrum = spectra[ibx]
        wn_spectrum = wn_spectra[ibx]

        # calculate the Fourier component for this baseline
        # (must be a better way of doing this)
        argx = numpy.radians(b_x * (wn * 100.0) * spatial_axis / 3600.0)
        argy = numpy.radians(b_y * (wn * 100.0) * spatial_axis / 3600.0)
        valx = numpy.exp(2.0j * math.pi * argx)
        valy = numpy.exp(2.0j * math.pi * argy)

        # get the visibility for this spectral plane
        vis = spectrum[wn_spectrum==wn]

        # calculate and coadd the contribution
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
        """Constructor.

        Parameters:
        previous_results - Current results structure of the simulation run.
        job_server       - ParallelPython job server.
        """
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()
        self.nuvspectra = None      

    def run(self):
#        print 'DirtyImage.run'

        # get relevant fts info
        fts = self.previous_results['fts']
        fts_wnmin = fts['wnmin']

        # get spatial info
        cubeparameters = self.previous_results['cubeparameters']
        npix = cubeparameters['npix']
        spatial_axis = cubeparameters['spatial axis [arcsec]']

        # get spectra at each baseline
        uvspectra = self.previous_results['reduceinterferogram']\
          ['baseline_uvspectrum'].values()
        self.nuvspectra = len(uvspectra)

        # get wavenumber axis
        wavenumber = uvspectra[0].wavenumber
        wavenumber = wavenumber[wavenumber>fts_wnmin]

        # dirty image cube
        dirtyimage = np.zeros([npix, npix, len(wavenumber)], np.float)
        dirtybeam = np.zeros([npix, npix, len(wavenumber)], np.float)

        # calculate dirty image for each wn
        # uvspectrum objects don't pickle which means they can't be
        # passed to the parallel processes as they are, so unpack them
        # first
        b_x_list = []
        b_y_list = []
        spectra = []
        wn_spectra = []
        for uvspectrum in uvspectra:
            b_x_list.append(uvspectrum.baseline_x)
            b_y_list.append(uvspectrum.baseline_y)
            spectra.append(uvspectrum.spectrum)
            wn_spectra.append(uvspectrum.wavenumber)

        # start a job for each wavenumber plane
        jobs = {}
        for iwn,wn in enumerate(wavenumber):
            # submit jobs
            indata = (b_x_list, b_y_list, spectra, wn_spectra, iwn, wn, 
              spatial_axis, npix,)
            jobs[wn] = self.job_server.submit(
              calculate_dirty_plane, 
              indata,
              (),
              ('numpy','math',))

        # collect and store results
        for iwn,wn in enumerate(wavenumber):
            if jobs[wn]() is None:
                raise Exception, 'calculate_dirty_plane has failed'

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
  cube        : {npix} x {npix} x {nwn}
  # baselines : {nbaselines}
'''.format(
           npix=np.shape(self.result['dirtyimage'])[0],
           nwn=np.shape(self.result['dirtyimage'])[2],
           nbaselines=self.nuvspectra)

