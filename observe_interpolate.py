from __future__ import absolute_import

import collections
import numpy as np
import numpy
import psutil
import time

import matplotlib.pyplot

def data_size(sky_cube, amplitude_beam_1, amplitude_beam_2):
    """Routine to calculate size of dominant data arrays in
    calculate_visibility.
    """
    result = sky_cube.nbytes + amplitude_beam_1.nbytes + \
      amplitude_beam_2.nbytes
    return result

def calculate_visibility(sky_cube, wn_axis, spatial_axis, 
  m1_area, td0L, td0R, amplitude_beam_1, amplitude_beam_2, obs_timeline,
  parallel):
    """Routine to calculate the visibility for a specified
    baseline for all planes in a sky model.

    Parameters:
    sky_cube          - 3d sky cube [i,j,wn]
    wn_axis           - wavenumbers of frequency axis in cm-1
    spatial_axis      - the spatial offsets of i,j axes
    m1_area           - area of collector primaries in m**2
    td0L              - transmission from primary to detector through
                      - left interferometer arm
    td0R              - transmission from primary to detector through
                      - right interferometer arm
    amplitude_beam_1  - the complex amplitude beam of collector 1
    amplitude_beam_2  - the complex amplitude beam of collector 2
    obs_timeline      - a dict containing the instrument configurations
                        to be simulated 
    parallel          - is this method being run in parallel with other
                        instances
    """ 
    nx,ny,nwn = numpy.shape(sky_cube)
    # spectra will hold the spectral points for each time and wn. The 
    # complete spectrum can be transformed (done outside this
    # routine) to give the interferogram at that point
    spectra = collections.defaultdict(dict)

    # i1 and i2 hold the DC component of the power falling onto the 
    # detector from each collector
    i1 = collections.defaultdict(dict) 
    i2 = collections.defaultdict(dict) 

    # width of spectral point in cm-1
    delta_wn = numpy.abs(wn_axis[1] - wn_axis[0])

    # find indices of centre of sky
    centre = numpy.argmin(numpy.abs(spatial_axis))

    # get grid of indices for sky, offset origin to centre
    grid = numpy.indices([nx, ny], numpy.float)
    jgrid = grid[0] - centre
    igrid = grid[1] - centre

    # calculate the u/v axis values (u axis == v axis in this case)
    u = numpy.fft.fftfreq(nx, spatial_axis[1] - spatial_axis[0])

    # initialize cache
    cache_bangle = None
    cache_bangle_limit = numpy.deg2rad(2.0)
    nmiss = 0
    ctime = 0.0

    obs_times = obs_timeline.keys()
    obs_times.sort()
    for it,t in enumerate(obs_times):
        if it%1000 == 0 and not parallel:
            print it
        config = obs_timeline[t]

        # if the baseline is bad then there is no signal
        if config.baseline_flag:
            for iwn, wn in enumerate(wn_axis):
                spectra[t][wn] = 0.0
            continue

        bx = config.baseline_x
        by = config.baseline_y
        bz = config.baseline_z

        bangle = numpy.arctan2(bx, by)
        if cache_bangle is None or \
          abs(cache_bangle - bangle) > cache_bangle_limit:
            cstart = time.clock()
            nmiss += 1
            #print 'calculating beam * sky because beam angle change has exceeded limit'
            #print wn_axis[0], t, bx, by, cache_bangle, bangle, cache_bangle_limit 
            # baseline angle has changed by more than specified limit,
            # recalculate the sky * conj(psf1) * psf2
        
            # ..calculate the jbeam,ibeam (=y,x) of the psf
            #   onto which each map j,i falls. 
            jbeam = jgrid * numpy.cos(bangle) + igrid * numpy.sin(bangle)
            ibeam = -jgrid * numpy.sin(bangle) + igrid * numpy.cos(bangle)  

            # ..remove origin offset
            jbeam += centre
            ibeam += centre
 
            # ..convert jbeam, ibeam to nearest integer values
            jbeam = numpy.rint(jbeam).astype(int)
            ibeam = numpy.rint(ibeam).astype(int)

            jbeam[jbeam<0] = 0
            jbeam[jbeam>nx-1] = nx-1
            ibeam[ibeam<0] = 0
            ibeam[ibeam>nx-1] = nx-1

            # following bit can be uncommented if you want to check
            # that the beam rotation code is working. It generates
            # plots that show the unrotated and rotated beams
            # together, with a horizontal line of zeros written
            # through the unrotated beam as a reference. 
            #res1 = numpy.abs(amplitude_beam_1[:,:,0])
            #res1[centre-1:centre+1,:] = 0.0
            #res = res1[jbeam, ibeam]
            #imshape = res.shape

            #plt = matplotlib.pyplot
            #plt.figure()

            #plt.subplot(211)
            #plt.imshow(res1, interpolation='nearest', origin='lower',
            #  aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
            #  vmax=numpy.max(res1[imshape[0]/4:imshape[0]*3/4,
            #  imshape[1]/4:imshape[1]*3/4]),
            #  vmin=numpy.min(res1[imshape[0]/4:imshape[0]*3/4,
            #  imshape[1]/4:imshape[1]*3/4]))
            #plt.colorbar(orientation='vertical')
            #plt.axis('image')
            #plt.title('basic')

            #plt.subplot(212)
            #plt.imshow(res, interpolation='nearest', origin='lower',
            #  aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
            #  vmax=numpy.max(res[imshape[0]/4:imshape[0]*3/4,
            #  imshape[1]/4:imshape[1]*3/4]),
            #  vmin=numpy.min(res[imshape[0]/4:imshape[0]*3/4,
            #  imshape[1]/4:imshape[1]*3/4]))
            #plt.colorbar(orientation='vertical')
            #plt.axis('image')
            #plt.title('bangle %s' % bangle)

            #filename = 'rottest%s%s.png' % (bangle, wn_axis[0])
            #plt.savefig(filename)
            #plt.close()

            # ..multiply the sky cube by the collectors' complex amp beams 
            # ..of each telescope
            # .. = Amp1 * td0L * Amp2 * td0R * delta_wn
            # ..where Amp1, Amp2 are sky cube amplitudes from telescopes 1 and 2
            # ..td0L and td0R are transmissions from telescope to detector 
            # ..L (telescope 1) and R (telescope 2) arms of interferometer.
            # ..Here, sky_cube is the sky intensity so 
            # .. Amp1 = sqrt (sky_cube1 * A1)
            # .. Amp2 = sqrt (sky_cube2 * A2)
            # ..where A1, A2 are the collector areas and for now 
            # ..sky_cube1 = sky_cube2 = sky_cube and A1 = A2 = ??. 
            f1 = sky_cube * m1_area * td0L * td0R * delta_wn * \
              numpy.conj(amplitude_beam_1[jbeam, ibeam]) * \
              amplitude_beam_2[jbeam, ibeam]
            sky_sum = numpy.sum(sky_cube, axis=(0,1))

            # ..calculate the FT of the sky planes
            # ....move centre of sky image to origin
            temp = numpy.fft.fftshift(f1, axes=(0,1,))
            # ....2d fft
            sky_fft = numpy.fft.fft2(temp, axes=(0,1,))

            cache_bangle = bangle
            ctime += (time.clock() - cstart)

        for iwn, wn in enumerate(wn_axis):
            lamb = 1.0 / (wn * 100.0)

            # find where the baseline falls on the sky_fft cube
            ang_freq = 1.0 / (numpy.rad2deg(lamb / bx) * 3600.0)
            col = ang_freq / (u[1] - u[0])

            ang_freq = 1.0 / (numpy.rad2deg(lamb / by) * 3600.0)
            row = ang_freq / (u[1] - u[0])

            col_lo = int(numpy.floor(col))
            col_hi = int(numpy.ceil(col))
            row_lo = int(numpy.floor(row))
            row_hi = int(numpy.ceil(row))

            # amps/phases at these angular freqs - using amp/phase
            # instead of real and imaginary components as these
            # should not suffer from wrapping issues
            amp_lolo = numpy.abs(sky_fft[row_lo,col_lo,iwn])
            amp_hilo = numpy.abs(sky_fft[row_hi,col_lo,iwn])
            amp_lohi = numpy.abs(sky_fft[row_lo,col_hi,iwn])
            amp_hihi = numpy.abs(sky_fft[row_hi,col_hi,iwn])

            ang_lolo = numpy.angle(sky_fft[row_lo,col_lo,iwn])
            ang_hilo = numpy.angle(sky_fft[row_hi,col_lo,iwn])
            ang_lohi = numpy.angle(sky_fft[row_lo,col_hi,iwn])
            ang_hihi = numpy.angle(sky_fft[row_hi,col_hi,iwn])

            # try to fix phase wrapping problems 
            unwrapped = numpy.unwrap(
              numpy.array([ang_lolo, ang_hilo, ang_lohi, ang_hihi]))
            ang_lolo = unwrapped[0] 
            ang_hilo = unwrapped[1]
            ang_lohi = unwrapped[2] 
            ang_hihi = unwrapped[3] 
 
            # interpolate the fft value using bilinear interpolation
            amp = amp_lolo + \
              (amp_lohi - amp_lolo) * (col - col_lo) + \
              (amp_hilo - amp_lolo) * (row - row_lo) + \
              (amp_lolo - amp_lohi - amp_hilo + amp_hihi) * \
              (col - col_lo) * (row - row_lo)

            angle = ang_lolo + \
              (ang_lohi - ang_lolo) * (col - col_lo) + \
              (ang_hilo - ang_lolo) * (row - row_lo) + \
              (ang_lolo - ang_lohi - ang_hilo + ang_hihi) * \
              (col - col_lo) * (row - row_lo)

            vis = amp * numpy.exp(1.0j * angle)
            spectra[t][wn] = vis
            i1[t][wn] = sky_sum[iwn] * m1_area * td0L * delta_wn
            i2[t][wn] = sky_sum[iwn] * m1_area * td0R * delta_wn

    print 'rotated beam recalculated %s times, time %s' % (nmiss, ctime)
    return spectra, i1, i2


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

    def run(self, parallel=True):
        """Method invoked to calculate the interferograms.
        """
        print 'Observe.run'
        print 'start', time.clock()

        # access telescope information
        telescope = self.previous_results['telescope']
        m1_diameter = telescope['m1_diameter']
        m1_area = np.pi * (m1_diameter/2)**2

        # access primary beam information
        primarybeams = self.previous_results['primarybeams']
        amp_beam_1 = primarybeams['collector 1 amplitude beam']
        amp_beam_2 = primarybeams['collector 2 amplitude beam']

        # access FTS information
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        # convert delta_opd from cm to m
        delta_opd = fts['delta_opd'] / 100.0
        smec_opd_to_mpd = fts['smec_opd_to_mpd']

        # access interferometer information
        noise = self.previous_results['noise']
        td0L = noise['td0L']
        td0R = noise['td0R']

        # access model sky
        skygenerator = self.previous_results['skymodel']
        sky_model = skygenerator['sky model']
        spatial_axis = self.result['spatial axis [arcsec]'] = \
          skygenerator['spatial axis [arcsec]']
        nx = len(spatial_axis)

        # get list of instrument configuratons       
        timeline = self.previous_results['timeline']
        obs_timeline = timeline['obs_timeline']
#        if timeline['pointing_error_type'].upper().strip() != 'ZERO':
#            print '..pointing errors are enabled, calculating the interferograms will take a long time'

        # calculate the measured result for each configuration
        previous_config = None
        observed_times = obs_timeline.keys()
        observed_times.sort()

        if not parallel:
            # do everything in one call to 'calculate visibility'
            chunks = []
            chunks.append(slice(0, len(fts_wn_truncated)))

        else:

            # decide how to slice up the problem frequency-wise
            ncpus = self.job_server.get_ncpus()        
            memory = psutil.virtual_memory()
            memory = memory.total
            dsize = data_size(sky_model, amp_beam_1.items()[0][1].data,
              amp_beam_2.items()[0][1].data)
            print 'ncpus=%s memory=%s data_size=%s' % (ncpus, memory, dsize)

            if dsize > memory:
                raise Exception, 'data size larger than physical memory, not handled'

            chunks = []
            nchunks = ncpus
            slice_size = len(fts_wn_truncated) / nchunks
    
            for chunk in range(nchunks):
                chunks.append(slice(chunk * slice_size, (chunk+1) * slice_size))

            # ensure last chunk covers required range
            last = chunks.pop()
            chunks.append(slice(last.start, len(fts_wn_truncated)))

        # now tackle the timeline
        # ..in preparation, get some info on the range of beam models available
        model_grid_1 = amp_beam_1.keys()
        baselines_1 = set()
        for b in model_grid_1:
            baselines_1.update([b])
        baselines_1 = list(baselines_1)
        baselines_1.sort()

        model_grid_2 = amp_beam_2.keys()
        baselines_2 = set()
        for b in model_grid_2:
            baselines_2.update([b])
        baselines_2 = list(baselines_2)
        baselines_2.sort()

        # work through the configurations to be observed
        istart = 0
        chunk_model_1 = None
        chunk_model_2 = None
        looping = True

        while looping:
        
            # assemble the times to go into the next timeline chunk
            timeline_chunk = []
            for ichunk in range(istart, len(observed_times)):

                # stop building if chunk len is 40000
                if len(timeline_chunk) > 40000:
                    istart = ichunk
                    break

                config = obs_timeline[observed_times[ichunk]]

                # ignore if baseline flagged bad
                if config.baseline_flag:
                    continue

                # stop building if either of the beam models change because 
                # of the varying baseline
                bx = config.baseline_x
                by = config.baseline_y
                bz = config.baseline_z
                blen = np.sqrt(bx**2 + by**2 +bz**2)

                model_b_1 = baselines_1[np.argmin(np.abs(baselines_1 - blen))]
                model_b_2 = baselines_2[np.argmin(np.abs(baselines_2 - blen))]

                if ichunk==istart:
                    chunk_model_1 = model_b_1
                    chunk_model_2 = model_b_2
                elif chunk_model_1 is not None and chunk_model_1 != model_b_1:
                    istart = ichunk
                    break
                elif chunk_model_2 is not None and chunk_model_2 != model_b_2:
                    istart = ichunk
                    break

                # OK, append this time to the chunk
                timeline_chunk.append(ichunk)
                looping = (ichunk != len(observed_times)-1)

            print '..timeline chunk', timeline_chunk[0], timeline_chunk[-1]

            obs_timeline_chunk = {}
            for it in timeline_chunk:
                t = observed_times[it]
                obs_timeline_chunk[t] = obs_timeline[t] 
                
            # get the beam models appropriate to this chunk
            amp_beam_1_chunk = amp_beam_1[model_b_1].data
            amp_beam_2_chunk = amp_beam_2[model_b_2].data
             
            if not parallel:
                chunk = chunks[0]
                job_id = (chunk.start, chunk.stop)
                powers = {}
                powers[job_id] = \
                  calculate_visibility(
                  sky_model[:,:,chunk],
                  fts_wn_truncated[chunk],
                  spatial_axis,
                  m1_area,
                  td0L,
                  td0R,
                  amp_beam_1_chunk[:,:,chunk], amp_beam_2_chunk[:,:,chunk], 
                  obs_timeline_chunk,
                  parallel=False)

            else:
                # submit jobs
                jobs = {}
                for chunk in chunks:
                    indata = (
                      sky_model[:,:,chunk],
                      fts_wn_truncated[chunk],
                      spatial_axis,
                      m1_area,
                      td0L,
                      td0R,
                      amp_beam_1_chunk[:,:,chunk],
                      amp_beam_2_chunk[:,:,chunk], 
                      obs_timeline_chunk,
                      True)

                    job_id = (chunk.start, chunk.stop)
                    print '....starting wavenumber chunk', job_id
                    jobs[job_id] = self.job_server.submit(calculate_visibility,
                      indata, (), ('numpy','collections','matplotlib.pyplot','time',))

                # collect results
                powers = {}
                i1 = {}
                i2 = {}
                for chunk in chunks:
                    job_id = (chunk.start, chunk.stop)
                    if jobs[job_id]() is None:
                        raise Exception, 'calculate_visibility has failed for planes %s' % str(job_id)
                    powers[job_id], i1[job_id], i2[job_id] = jobs[job_id]()

            # assemble spectra from single/parallel results
            keys = powers.keys()
            spectra = powers[keys[0]]
            i1_spectra = i1[keys[0]]
            i2_spectra = i2[keys[0]]
            for k in keys[1:]:
                for t in spectra.keys():
                    spectra[t].update(powers[k][t])
                    i1_spectra[t].update(i1[k][t])
                    i2_spectra[t].update(i2[k][t])

            # calculate power for each time
            for t in obs_timeline_chunk:
                config = obs_timeline[t]
                opd = config.smec_position / smec_opd_to_mpd
                opd_ipos = opd / delta_opd

                spectrum = np.zeros(np.shape(fts_wn), np.complex)
                for iwn, wn in enumerate(fts_wn):
                    try:
                        spectrum[iwn] = spectra[t][wn]
                    except:
                        spectrum[iwn] = 0.0

                # I'll try to add some theoretical background here,
                # for my own benefit on future visits if nothing else.
                # Neglecting transmission losses, at the detector:
                # 
                #  amp = A1.exp(i.omega.t) + A2.exp(i.([b.theta+l/lambda] + omega).t)
                #
                # where b.theta is the delay due baseline projection, 
                # l is the opd due to the FTS mirror.
                #
                #    I = amp.amp*
                # 
                # where * denotes complex conjugate. In principle, A1 and A2 are 
                # complex and functions of u,v and lambda.
                # This gives:  
                #    I = I1 + I2 + A1.A2*.exp(-i.b.theta/lambda).exp(-i.l/lambda) +
                #                  A1*.A2.exp(i.b.theta/lambda).exp(i.l/lambda)
                #
                # I think the last 2 terms if integrated over + and - frequency
                # axes (frequency = c/lambda) are the Fourier Transform of
                # function that is A1.A2*.exp(-i.b.theta/lambda) on +ve axis,
                # A1*.A2.exp(i.b.theta/lambda) on -ve axis. The function being
                # transformed is Hermitian, hence the result is real.
                #
                # If A1=A2 then the exponentials simplify to 2*cosine. However,
                # the method used here is more general.
                # (I hope this is right - doubts about that are one reason for this
                # comment).
                # 
                # Additional complication. The mesh filter used as a beamsplitter
                # in the FTS adds a pi/2 phase shift to the reflected beam. Thus
                # we have:
                #
                #  amp = A1.exp(i.omega.t) + i.A2.exp(i.([b.theta+l/lambda] + omega).t)
                #
                # at the start of the theory above. This means that the 
                # Hermitian function is i.A1.A2*.exp(-i.b.theta/lambda) on +ve axis,
                # -i.A1*.A2.exp(i.b.theta/lambda) on -ve axis.

                # calculate shift needed to move point at opd to 0
                nfreq = len(fts_wn)
                ndoublefreq = 2 * (nfreq - 1)
                shift = np.zeros([ndoublefreq], dtype=np.complex)
                shift[:nfreq] = np.arange(nfreq, dtype=np.complex)
                shift[nfreq:] = -np.arange(nfreq, dtype=np.complex)[-2:0:-1]
                shift = np.exp((2.0j * np.pi * opd_ipos * shift) / \
                  float(ndoublefreq))

                # multiply by i because of phase shift at metal-mesh
                # beamsplitter 
                spectrum *= 1j
                # reflect conjugate spectrum about 0 to give Hermitian version
                reflected_spectrum = np.zeros([ndoublefreq], np.complex)
                reflected_spectrum[:nfreq] = spectrum
                reflected_spectrum[nfreq:] = np.conjugate(spectrum[-2:0:-1])

                # apply phase shift and calculate inverse fft
                reflected_spectrum *= shift

                spectrum_fft = np.fft.ifft(reflected_spectrum)
                # remove the 1/n normalisation of ifft where n is the
                # length of spectrum not reflected_spectrum
                spectrum_fft *= ndoublefreq

                # interferogram = I1 + I2 + FT(Hermitian spectrum)
                # I1 and I2 can vary with t if the beam of either collector
                # moves about
                power = np.sum(i1_spectra[t].values()) + \
                  np.sum(i2_spectra[t].values()) + \
                  spectrum_fft[0].real

                # set the measured value
                config = config._replace(data = power, pure_data=power)
                obs_timeline[t] = config

        self.result['observed_timeline'] = obs_timeline
        print 'stop', time.clock()
        return self.result

    def __repr__(self):
        return '''
Observe:
  timelength length: {timeline_len}
'''.format(
          timeline_len=len(self.result['observed_timeline']))

