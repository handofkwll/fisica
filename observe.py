from __future__ import absolute_import

import collections
import numpy as np
import numpy
import psutil
import time
import visibility

#used for debugging
#import matplotlib.pyplot as plt
#import common.commonobjects as co

def data_size(sky_cube, amplitude_beam_1, amplitude_beam_2):
    """Routine to calculate size of dominant data arrays in
    calculate_visibility.
    """
    result = sky_cube.nbytes + amplitude_beam_1.nbytes + \
      amplitude_beam_2.nbytes
    return result

def calculate_visibility(smec_opd_to_mpd, sky_cube, wn_axis,
  spatial_axis,
  amplitude_beam_1, amplitude_beam_2, obs_timeline, times):
    """Routine to calculate the visibility for a specified
    baseline for one plane in a sky model.

    Parameters:
    baseline          - (u,v) in metres
    wn                - wavenumber of observation (cm-1)
    sky_plane         - 2d sky image at this wn
    spatial_freq_axis - The spatial frequency axis of the sky plane
                        Fourier transform
    """ 
    def rounded(in_tuple):
        """Routine to round the elements of a tuple being used
        as a dict key. The idea is to ensure that keys that
        are close to identical are treated as identical.
        """
        out = [int(round(1000000 * v)) for v in in_tuple]
        return tuple(out) 

    nx,ny,nwn = numpy.shape(sky_cube)
    power = {}
    angles = numpy.radians(spatial_axis / 3600.0) 

    dpower_cache = {}
    ncached = 0
    ntotal = 0

    f1_cache = {}
    nsearched_f1 = 0
    nhit_f1 = 0

    f1timesf2_cache = {}
    nsearched_baseline = 0
    nhit_baseline = 0

    # using Anthony Murphy's formalism throughout

    # first, calculate f1 = B * conj(psf1) * psf2
    # only cache 0 pointing error case
    # TBD - shift psfs according to pointing errors   
    psf1 = amplitude_beam_1
    psf2 = amplitude_beam_2

#    f1_cache[rounded((0.0,0.0,0.0,0.0))] = \
#      sky_cube * numpy.conj(psf1) * psf2

    for it,t in enumerate(times):
        #print it
        config = obs_timeline[t]
        opd = config.smec_position / smec_opd_to_mpd
        ntotal += 1
       
        power[t] = 0.0

        # first, calculate f1 = B * conj(psf1) * psf2
#        nsearched_f1 += 1
#        f1_key = rounded((
#          config.pointing1_x,
#          config.pointing1_y,
#          config.pointing2_x,
#          config.pointing2_y))
#        f1_cube = f1_cache.get(f1_key)
#        if f1_cube is None:
#            f1_cube = sky_cube * numpy.conj(psf1) * psf2
#        else:
#            nhit_f1 += 1

        for iwn, wn in enumerate(wn_axis):

            # see if the result is in the cache
      
            if dpower_cache.has_key(
              rounded((opd, wn, 
              config.pointing1_x, config.pointing1_y,
              config.pointing2_x, config.pointing2_y,
              config.baseline_x, config.baseline_y))):

                dpower = dpower_cache[rounded((opd, wn, 
                  config.pointing1_x, config.pointing1_y,
                  config.pointing2_x, config.pointing2_y,
                  config.baseline_x, config.baseline_y))]

                power[t] += dpower.real
                # purpose of this commented out section
                # described below
#                if iwn==len(wn_axis)-1:
#                    power[t] += dpower.real
#                else:
#                    power[t] += 2.0 * dpower.real

                ncached += 1
                continue

            # see if the required f1timesf2 is in cache

            f1timesf2_key = rounded((wn,
              config.pointing1_x, config.pointing1_y,
              config.pointing2_x, config.pointing2_y,
              config.baseline_x, config.baseline_y))
            
            f1timesf2 = f1timesf2_cache.get(f1timesf2_key)
            if f1timesf2 is None:

                # calculate it

                # see if cached f1 = B * conj(psf1) * psf2
                nsearched_f1 += 1
                f1_key = rounded((wn,
                  config.pointing1_x,
                  config.pointing1_y,
                  config.pointing2_x,
                  config.pointing2_y))
                f1_plane = f1_cache.get(f1_key)
                if f1_plane is None:
                    f1_plane = sky_cube[:,:,iwn] * \
                      numpy.conj(psf1[:,:,iwn]) * psf2[:,:,iwn]
                else:
                    nhit_f1 += 1

                #  calculate f2 = exp(-j * k(theta,phi).b)
                #  only cache range of wn for one baseline at a time
                f2x = numpy.exp(2.0j * numpy.pi * wn * 100.0 * (-angles * config.baseline_x))
                f2y = numpy.exp(2.0j * numpy.pi * wn * 100.0 * (-angles * config.baseline_y))
                f2 = numpy.ones([nx,nx], numpy.complex)
                # using array broadcasting to use C indexing in numpy and not
                # in Python
                f2 *= f2x
                f2shape = list(numpy.shape(f2))
                f2shape.reverse()
                f2.reshape(f2shape)
                f2 *= f2y
                f2shape.reverse()
                f2.reshape(f2shape)

                #  calculate 2 * sum[B * exp(-j * k(theta,phi).b)]
                f1timesf2 = 2.0 * numpy.sum(f2 * f1_plane)
                f1timesf2_cache[rounded((wn,
                  config.pointing1_x, config.pointing1_y,
                  config.pointing2_x, config.pointing2_y,
                  config.baseline_x, config.baseline_y))] = f1timesf2

            #  calculate 2 * sum[B * exp(j * chi(theta, phi, tau)) * conj(psf1) * psf2]
            delta_power = f1timesf2 * numpy.exp(2.0j * numpy.pi * wn * 100.0 * opd)

            power[t] += delta_power.real
            # following section to handle FFT problem (not
            # sure its right to call it aliasing) that occurs if there
            # is flux at the maximum frequency. 
            #
            # The interferogram calculated here is not the same
            # as that via numpy FFT unless there is zero flux at
            # maximum frequency. Numpy FFT
            # for even n effectively counts each frequency twice
            # but the highest frequency once (for real data).
            #
            # For now enforce zero flux at wnmax by making
            # cutoffmax less than wnmax.
            # - at least, that is what I think is happening
#            if iwn==len(wn_axis)-1:
#            else:
#                power[t] += 2.0 * delta_power.real

            dpower_cache[rounded((opd, wn,
              config.pointing1_x, config.pointing1_y,
              config.pointing2_x, config.pointing2_y,
              config.baseline_x, config.baseline_y))] = delta_power.real

            #print wn, f1timesf2.real

#           debug plotting
#            plt.figure()
#            plt.imshow(factor.real, interpolation='nearest', origin='lower',
#              aspect='equal')
#            plt.colorbar(orientation='vertical')
#            plt.axis('image')
#            filename = 'debug.png'
#            plt.savefig(filename)
#            plt.close()
#            x = 1/0

    fraction_cached = float(ncached)/float(ntotal)

#    # data are pickled for transfer between ParallelPython processes
#    # and this won't work for defaultdict objects - so convert to dict
#    return dict(power), fraction_cached
    f1_hit_rate = float(nhit_f1) / float(nsearched_f1)
    return power, fraction_cached, f1_hit_rate

def fftshift(data, shift):
    """Method to shift a 2d complex data array by applying the
       given phase shift to its Fourier transform.

    Parameters:
    data  - 2d data to be shifted
    shift - 2d array containing phase shift
    """
    # move centre of image to array origin
    temp = numpy.fft.fftshift(data)
    # 2d fft
    temp = numpy.fft.fft2(temp)
    # apply phase shift
    temp *= shift
    # transform and shift back
    temp = numpy.fft.ifft2(temp)
    temp = numpy.fft.fftshift(temp)

    return temp

def fields_match(previous_config, config, exclude_fields=[]):
    """Function to compare 2 configs, ignoring any differences in
    the components listed in 'exclude_fields'.

    Parameters:
    previous_config - namedtuple containing the previous configuration
    config          - namedtuple with the configuration to compare
    exclude_fields  - list of tuple fields to be excluded in the
                      comparison
    """
    if previous_config is None:
        return False

    if not exclude_fields:
        return previous_config==config
    else:
        # convert to dicts and eliminate exclude_fields
        previous_dict = previous_config._asdict()
        current_dict = config._asdict()
        for field in exclude_fields:
            ignore = previous_dict.pop(field, None)
            ignore = current_dict.pop(field, None)

        previous_set = set([k for k in previous_dict.items()])
        current_set = set([k for k in current_dict.items()])
        diff = previous_set.symmetric_difference(current_set)

        debug = False
        if(debug):
            print diff

        return not(bool(diff))


class Observe(object):
    """Class to compute interferograms.
    """

    def __init__(self, parameters, previous_results, job_server):
        """Constructor.

        Parameters:
        parameters       - Dict with parameters from the Excel files.
        previous_results - Current results structure of the simulation run.
        job_server       - ParallelPython job server. 
        """
        self.parameters = parameters
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        """Method invoked to calculate the interferograms.
        """
        print 'Observe.run'
        print time.clock()

        # access primary beam information
        primarybeams = self.previous_results['primarybeams']
        amp_beam_1 = primarybeams['primary amplitude beam'].data
        amp_beam_2 = primarybeams['primary amplitude beam'].data

        # access required FTS information
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        # convert delta_opd from cm to m
        delta_opd = fts['delta_opd'] / 100.0
        smec_opd_to_mpd = fts['smec_opd_to_mpd']

        # access to model sky
        skygenerator = self.previous_results['skygenerator']
        sky_model = skygenerator['sky model']
        spatial_axis = self.result['spatial axis'] = skygenerator['spatial axis']
        nx = len(spatial_axis)

        # get list of instrument configuratons       
        timeline = self.previous_results['timeline']
        obs_timeline = timeline['obs_timeline']

        # calculate the measured result for each configuration
        previous_config = None
        observed_times = obs_timeline.keys()
        observed_times.sort()

        debug_plotted = False

        # decide how to slice up the problem
        ncpus = self.job_server.get_ncpus()        
        memory = psutil.virtual_memory()
        memory = memory.total
        dsize = data_size(sky_model, amp_beam_1, amp_beam_2)
        print 'ncpus=%s memory=%s data_size=%s' % (ncpus, memory, dsize)

        if dsize > memory:
            raise Exception, 'data size larger than physical memory, not handled'

        chunks = []
        nchunks = ncpus
        slice_size = len(fts_wn_truncated) / nchunks

        for chunk in range(nchunks):
            chunks.append(slice(chunk * slice_size, (chunk+1) * slice_size))

        last = chunks.pop()
        chunks.append(slice(last.start, len(fts_wn_truncated)))

        # direct call that can be used for debugging
        chunk = 0
        chunks[0] = slice(0, len(fts_wn_truncated))
        job_id = (chunks[chunk].start, chunks[chunk].stop)
        powers = {}
        fraction_cached = {}
        f1_hit_rate = {}
        powers[job_id], fraction_cached[job_id], f1_hit_rate[job_id] = \
          calculate_visibility(
          smec_opd_to_mpd,
          sky_model[:,:,chunks[chunk]],
          fts_wn_truncated[chunks[chunk]], spatial_axis,
          amp_beam_1[:,:,chunks[chunk]], amp_beam_2[:,:,chunks[chunk]], 
          obs_timeline, observed_times)

        # submit jobs
#        jobs = {}
#        for chunk in chunks[:4]:
#            indata = (smec_opd_to_mpd,
#                      sky_model[:,:,chunk],
#                      fts_wn_truncated[chunk],
#                      spatial_axis,
#                      amp_beam_1[:,:,chunk],
#                      amp_beam_2[:,:,chunk], 
#                      obs_timeline, observed_times,)

#            job_id = (chunk.start, chunk.stop)
#            print 'starting ', job_id
#            jobs[job_id] = self.job_server.submit(calculate_visibility,
#              indata, (), ('numpy','collections',))

        # collect results
#        powers = {}
#        fraction_cached = {}
#        f1_hit_rate = {}
#        for chunk in chunks:
#            job_id = (chunk.start, chunk.stop)
#            if jobs[job_id]() is None:
#                raise Exception, 'calculate_visibility has failed for planes %s' % str(job_id)
#            powers[job_id], fraction_cached[job_id], f1_hit_rate[job_id] = jobs[job_id]()

        keys = powers.keys()
        visibilities = powers[keys[0]]
        print keys[0], 'fraction cached', fraction_cached[keys[0]], \
          f1_hit_rate[keys[0]]
        for k in keys[1:]:
            for t in visibilities.keys():
                visibilities[t] += powers[k][t]
#                print k, 'fraction cached', fraction_cached[k], f1_hit_rate[k]

        for t in observed_times:
            config = obs_timeline[t]
            config = config._replace(data = visibilities[t])
            obs_timeline[t] = config

        self.result['observed_timeline'] = obs_timeline
        print 'stop'
        print time.clock()
        return self.result

    def __repr__(self):
        return '''
Observe:
  timelength length: {timeline_len}
'''.format(
          timeline_len=len(self.result['observed_timeline']))

