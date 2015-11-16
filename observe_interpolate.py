from __future__ import absolute_import

import collections
import numpy as np
import numpy
import psutil
import time

#import matplotlib.pyplot
#import scipy.interpolate

def data_size(sky_cube, amplitude_beam_1, amplitude_beam_2):
    """Routine to calculate size of dominant data arrays in
    calculate_visibility.
    """
    # multiplier must match the amount that the sky cube is 
    # padded out by for the fft = calculate_visibility.pad_factor**2
    result = 16 * sky_cube.nbytes + amplitude_beam_1.nbytes + \
      amplitude_beam_2.nbytes
    return result

def calculate_visibility(sky_cube, wn_axis, spatial_axis, 
  m1_area, td0L, td0R, amplitude_beam_1, amplitude_beam_2, 
  beam_angle, pointing_offset_1, pointing_offset_2, obs_timeline,
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
    beam_angle        - the angle to which the beams are to be rotated
    pointing_offset_1 - tuple containing pointing offset in (x,y) for
                        collector 1
    pointing_offset_2 - tuple containing pointing offset in (x,y) for
                        collector 2
    obs_timeline      - a dict containing the instrument configurations
                        to be simulated 
    parallel          - is this method being run in parallel with other
                        instances
    """ 
    nx,ny,nwn = numpy.shape(sky_cube)

    # spectra will hold the spectral points for each time and wn. The 
    # complex spectrum can be transformed (done outside this
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
    grid = numpy.indices([ny, nx], numpy.float)
    jgrid = grid[0] - centre
    igrid = grid[1] - centre

    # calculate the u/v axis values (u axis == v axis in this case)
    pad_factor = 4
    u = numpy.fft.fftfreq(pad_factor * nx, 
      spatial_axis[1] - spatial_axis[0])
    u = numpy.fft.fftshift(u)

    obs_times = obs_timeline.keys()
    obs_times.sort()
    first = True

    for it,t in enumerate(obs_times):
        if it%1000 == 0 and not parallel:
            print it
        config = obs_timeline[t]

        # if the baseline is bad then there is no signal
        if config.baseline_flag:
            for wn in wn_axis:
                spectra[t][wn] = 0.0
                i1[t][wn] = 0.0
                i2[t][wn] = 0.0
            continue

        bx = config.baseline_x
        by = config.baseline_y
        bz = config.baseline_z

        if first:
            first = False

            # calculate the sky * conj(psf1) * psf2
            # ..calculate the jbeam,ibeam (=y,x) of the psf
            # ..onto which each map j,i falls, taking into account
            # ..both the beam rotation and the pointing offset
            bangle = numpy.deg2rad(beam_angle) 
            jbeam = jgrid * numpy.cos(bangle) + igrid * numpy.sin(bangle)
            ibeam = -jgrid * numpy.sin(bangle) + igrid * numpy.cos(bangle)  

            # ..pointing offset for collector 1
            jpointing1 = pointing_offset_1[0] / (spatial_axis[1] - spatial_axis[0])
            ipointing1 = pointing_offset_1[1] / (spatial_axis[1] - spatial_axis[0])
            jbeam1 = jbeam - (jpointing1 * numpy.cos(bangle) + ipointing1 * numpy.sin(bangle))
            ibeam1 = ibeam - (-jpointing1 * numpy.sin(bangle) + ipointing1 * numpy.cos(bangle))  

            # ..remove origin offset
            jbeam1 += centre
            ibeam1 += centre
 
            # ..convert jbeam, ibeam to nearest integer values
            jbeam1 = numpy.rint(jbeam1).astype(int)
            ibeam1 = numpy.rint(ibeam1).astype(int)

            jbeam1[jbeam1<0] = 0
            jbeam1[jbeam1>nx-1] = nx-1
            ibeam1[ibeam1<0] = 0
            ibeam1[ibeam1>nx-1] = nx-1

            # same for collector 2
            jpointing2 = pointing_offset_2[0] / (spatial_axis[1] - spatial_axis[0])
            ipointing2 = pointing_offset_2[1] / (spatial_axis[1] - spatial_axis[0])
            jbeam2 = jbeam - (jpointing2 * numpy.cos(bangle) + ipointing2 * numpy.sin(bangle))
            ibeam2 = ibeam - (-jpointing2 * numpy.sin(bangle) + ipointing2 * numpy.cos(bangle))  
            jbeam2 += centre
            ibeam2 += centre
            jbeam2 = numpy.rint(jbeam2).astype(int)
            ibeam2 = numpy.rint(ibeam2).astype(int)
            jbeam2[jbeam2<0] = 0
            jbeam2[jbeam2>nx-1] = nx-1
            ibeam2[ibeam2<0] = 0
            ibeam2[ibeam2>nx-1] = nx-1

            # following bit can be uncommented if you want to check
            # that the beam rotation code is working. It generates
            # plots that show the unrotated and rotated beams
            # together, with a horizontal line of zeros written
            # through the unrotated beam as a reference. 
            #res1 = numpy.abs(amplitude_beam_1[:,:,0])
            #res1[centre-1:centre+1,:] = 0.0
            #res1[centre-1:centre+1,centre-1:centre+1] = 2.0
            #res = res1[jbeam1, ibeam1]
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
            #plt.title('bangle %s pointing j %s i %s' % (bangle, jpointing, ipointing))

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
              numpy.conj(amplitude_beam_1[jbeam1, ibeam1]) * \
              amplitude_beam_2[jbeam2, ibeam2]
            sky_sum_1 = numpy.sum(
              sky_cube * numpy.abs(amplitude_beam_1[jbeam1,ibeam1]), axis=(0,1))
            sky_sum_2 = numpy.sum(
              sky_cube * numpy.abs(amplitude_beam_2[jbeam2,ibeam2]), axis=(0,1))

            # ..calculate the FT of the sky planes
            # ....move centre of sky image to origin
            # note to self, axis 2 is fastest changing in this array,
            # the ordering is not optimal
            padding = (pad_factor-1) * nx / 2
            f1 = numpy.pad(f1, ((padding, padding), (padding, padding), (0,0)),
              'constant', constant_values=0) 
            f1 = numpy.fft.fftshift(f1, axes=(0,1,))
            # ....2d fft
            sky_fft = numpy.fft.fft2(f1, axes=(0,1,))
            del f1
            sky_fft = numpy.fft.fftshift(sky_fft, axes=(0,1,))

            interp = {}
            for iwn, wn in enumerate(wn_axis):
                interp[iwn] = (
                  scipy.interpolate.RectBivariateSpline(u, u, sky_fft[:,:,iwn].real),        
                  scipy.interpolate.RectBivariateSpline(u, u, sky_fft[:,:,iwn].imag))        

        ang_freq_u = numpy.zeros(wn_axis.shape)
        ang_freq_v = numpy.zeros(wn_axis.shape)
        for iwn, wn in enumerate(wn_axis):
            # find where the baseline falls on the sky_fft cube
            lamb = 1.0 / (wn * 100.0)
            ang_freq_u[iwn] = 1.0 / (numpy.rad2deg(lamb / bx) * 3600.0)
            ang_freq_v[iwn] = 1.0 / (numpy.rad2deg(lamb / by) * 3600.0)

        vis_real = interp[iwn][0](ang_freq_v, ang_freq_u, grid=False)
        vis_imag = interp[iwn][1](ang_freq_v, ang_freq_u, grid=False)
        vis = vis_real + 1j * vis_imag

        for iwn, wn in enumerate(wn_axis):
            spectra[t][wn] = vis[iwn]
            i1[t][wn] = sky_sum_1[iwn] * m1_area * td0L * delta_wn
            i2[t][wn] = sky_sum_2[iwn] * m1_area * td0R * delta_wn

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
        beam_angle_res = min(
          primarybeams['collector 1 rotation resolution'],
          primarybeams['collector 2 rotation resolution'])

        # access FTS information
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        # convert delta_opd from cm to m
        delta_opd = fts['delta_opd'] / 100.0
        smec_opd_to_mpd = fts['smec_opd_to_mpd']

        # access interferometer information
        backgroundnoise = self.previous_results['backgroundnoise']
        td0L = backgroundnoise['td0L']
        td0R = backgroundnoise['td0R']

        # access model sky
        skygenerator = self.previous_results['skymodel']
        sky_model = skygenerator['sky model']
        spatial_axis = self.result['spatial axis [arcsec]'] = \
          skygenerator['spatial axis [arcsec]']
        nx = len(spatial_axis)

        # get list of instrument configuratons       
        timeline = self.previous_results['timeline']
        obs_timeline = timeline['obs_timeline']

        # calculate the measured result for each configuration
        previous_config = None
        observed_times = obs_timeline.keys()
        observed_times.sort()

        # decide how to slice up the problem frequency-wise
        ncpus = self.job_server.get_ncpus()        
        memory = psutil.virtual_memory()
        memory = memory.total
        dsize = data_size(sky_model, amp_beam_1.items()[0][1].data,
          amp_beam_2.items()[0][1].data)
        print 'ncpus=%s memory=%s data_size=%s' % (ncpus, memory, dsize)

        if dsize > memory:
            raise Exception, 'data size larger than physical memory, not handled'

        # ideally ncpus * memory per chunk ~ 0.1 * memory
        size_per_cpu = 0.125 * memory / ncpus
        #size_per_cpu = 0.25 * memory / ncpus
        nchunks = max(ncpus, int(np.ceil(dsize / size_per_cpu)))

        chunks = []
        slice_size = len(fts_wn_truncated) / nchunks
    
        for chunk in range(nchunks):
            chunks.append(slice(chunk * slice_size, (chunk+1) * slice_size))

        # ensure last chunk covers required range
        last = chunks.pop()
        chunks.append(slice(last.start, len(fts_wn_truncated)))

        print '..calculation will break frequency band into', nchunks, 'chunks'

        # decide how to slice up the problem beam-rotation wise
        # ..break 360 degrees up into beam angular resolution pieces
        print '..calculation will rotate primary beam at', beam_angle_res,\
          'deg intervals'
        n_beam_angles = int(360/beam_angle_res)
        exact_beam_angle_res = 360.0 / n_beam_angles
        beam_angles = np.arange(n_beam_angles) * exact_beam_angle_res
        beam_bins = {}
        for beam_angle in beam_angles:
            beam_bins[beam_angle] = []

        # ..fill segments with times where the segment is occupied
        for t in observed_times:
            config = obs_timeline[t]
            bx = config.baseline_x
            by = config.baseline_y
            bz = config.baseline_z

            bangle = numpy.arctan2(bx, by) * 180.0 / np.pi
            # each bin centre is effectively (bin_angle + 0.5 * exact_beam_angle_res)
            bangle_index = int(np.floor(bangle / exact_beam_angle_res))
            if bangle_index < 0:
                bangle_index += n_beam_angles
            bin_angle = bangle_index * exact_beam_angle_res
            beam_bins[bin_angle].append(t)

        #print 'beam_bins'
        #for k,v in beam_bins.iteritems():
        #    print k, len(v)

        # slice the problem up pointing wise
        # ..break pointing into 7x7 grid
        pointing_bins = collections.defaultdict(list)
        max_point = 0.0
        for t in observed_times:
            config = obs_timeline[t]
            pointing1_x = config.pointing1_x
            pointing1_y = config.pointing1_y

            max_point = max(max_point, pointing1_x, pointing1_y)

        # most pointings within 3 arcsec radius of nominal. Outer
        # ring has radius 4 * pfactor = 0.3
        pfactor = 0.3 / 4
        first_circle = [(2*pfactor, 0),
                        (2*pfactor*np.cos(45), 2*pfactor*np.sin(45)), 
                        (0, 2*pfactor), 
                        (2*pfactor*np.cos(135), 2*pfactor*np.sin(135)), 
                        (-2*pfactor, 0),
                        (2*pfactor*np.cos(-135), 2*pfactor*np.sin(-135)),
                        (0, -2*pfactor),
                        (2*pfactor*np.cos(-45), 2*pfactor*np.sin(-45))] 
        second_circle = [(4*pfactor, 0),
                         (4*pfactor*np.cos(22.5), 4*pfactor*np.sin(22.5)), 
                         (4*pfactor*np.cos(45), 4*pfactor*np.sin(45)), 
                         (4*pfactor*np.cos(67.5), 4*pfactor*np.sin(67.5)), 
                         (0, 4*pfactor), 
                         (4*pfactor*np.cos(112.5), 4*pfactor*np.sin(112.5)), 
                         (4*pfactor*np.cos(135), 4*pfactor*np.sin(135)), 
                         (4*pfactor*np.cos(157.5), 4*pfactor*np.sin(157.5)), 
                         (-4*pfactor, 0),
                         (4*pfactor*np.cos(-157.5), 4*pfactor*np.sin(-157.5)), 
                         (4*pfactor*np.cos(-135), 4*pfactor*np.sin(-135)),
                         (4*pfactor*np.cos(-112.5), 4*pfactor*np.sin(-112.5)), 
                         (0, -4*pfactor),
                         (4*pfactor*np.cos(-67.5), 4*pfactor*np.sin(-67.5)), 
                         (4*pfactor*np.cos(-45), 4*pfactor*np.sin(-45)), 
                         (4*pfactor*np.cos(-22.5), 4*pfactor*np.sin(-22.5))]
 
        for t in observed_times:
            config = obs_timeline[t]
            pointing1_x = config.pointing1_x
            pointing1_y = config.pointing1_y

            pointing_r = np.sqrt(pointing1_x**2 + pointing1_y**2)
            if pointing_r < pfactor:
                pointing_bins[0,0].append(t)
            elif pointing_r < 3*pfactor:
                pointing_ang = np.arctan2(pointing1_x, pointing1_y) * 180.0 / np.pi
                pointing_ang_ind = int(np.floor(pointing_ang / 45.0))
                if pointing_ang_ind < 0:
                    pointing_ang_ind += 8
                pointing_bins[first_circle[pointing_ang_ind]].append(t)
            else:
                pointing_ang = np.arctan2(pointing1_x, pointing1_y) * 180.0 / np.pi
                pointing_ang_ind = int(np.floor(pointing_ang / 22.5))
                if pointing_ang_ind < 0:
                    pointing_ang_ind += 16
                pointing_bins[second_circle[pointing_ang_ind]].append(t)
            
        # slice up the problem baseline-length-wise
        #pointing_bins[(0,0)] = observed_times

        # ..populate pointing bins with times 

        # assemble a list of wn_chunk, beam bin, pointing bin
        calculation_list = []
        for beam_angle,beam_bin in beam_bins.items():
            for pointing_offset_1,pointing_bin in pointing_bins.items():
                times = set(beam_bin)
                times = times.intersection(pointing_bin)
                times = list(times)
                if times:
                    # break times up into shorter chunks in crude
                    # avoid memory problems
                    time_chunks = [times[i:i+100000] for i in range(0, len(times), 100000)]
                    for time_chunk in time_chunks:
                        band_chunks = []
                        for chunk in chunks:
                            band_chunks.append((time_chunk,
                              beam_angle + 0.5 * exact_beam_angle_res,
                              pointing_offset_1,
                              chunk))
                        calculation_list.append(band_chunks)

        print '..calculation broken into', len(calculation_list), 'pieces' 

        # now tackle the timeline
        # ..in preparation, get some info on the range of beam models available
        model_grid_1 = amp_beam_1.keys()
#        baselines_1 = set()
#        for b in model_grid_1:
#            baselines_1.update([b])
#        baselines_1 = list(baselines_1)
#        baselines_1.sort()

#        model_grid_2 = amp_beam_2.keys()
#        baselines_2 = set()
#        for b in model_grid_2:
#            baselines_2.update([b])
#        baselines_2 = list(baselines_2)
#        baselines_2.sort()

        # work through the configurations to be observed
#        istart = 0
#        chunk_model_1 = None
#        chunk_model_2 = None
#        looping = True

#        while looping:
        
            # assemble the times to go into the next timeline chunk
#            timeline_chunk = []
#            for ichunk in range(istart, len(observed_times)):

                # stop building if chunk len is 40000
#                if len(timeline_chunk) > 40000:
#                    istart = ichunk
#                    break

#                config = obs_timeline[observed_times[ichunk]]

                # ignore if baseline flagged bad
#                if config.baseline_flag:
#                    continue

                # stop building if either of the beam models change because 
                # of the varying baseline
#                bx = config.baseline_x
#                by = config.baseline_y
#                bz = config.baseline_z
#                blen = np.sqrt(bx**2 + by**2 +bz**2)

#                model_b_1 = baselines_1[np.argmin(np.abs(baselines_1 - blen))]
#                model_b_2 = baselines_2[np.argmin(np.abs(baselines_2 - blen))]

#                if ichunk==istart:
#                    chunk_model_1 = model_b_1
#                    chunk_model_2 = model_b_2
#                elif chunk_model_1 is not None and chunk_model_1 != model_b_1:
#                    istart = ichunk
#                    break
#                elif chunk_model_2 is not None and chunk_model_2 != model_b_2:
#                    istart = ichunk
#                    break

                # OK, append this time to the chunk
#                timeline_chunk.append(ichunk)
#                looping = (ichunk != len(observed_times)-1)

#            print '..timeline chunk', timeline_chunk[0], timeline_chunk[-1]

        # submit jobs
        jobs = collections.deque()
        for counter,band_chunks in enumerate(calculation_list):
            print '..calculation piece', counter
            spectra = collections.defaultdict(dict)
            i1_spectra = collections.defaultdict(dict)
            i2_spectra = collections.defaultdict(dict)

            for calculation in band_chunks:
                times = calculation[0]
                beam_angle = calculation[1]
                pointing_offset_1 = calculation[2]
                wn_chunk = calculation[3]

                obs_timeline_chunk = {}
                for t in times:
                    obs_timeline_chunk[t] = obs_timeline[t] 

                # get the beam models appropriate to this chunk
                model_b_1 = amp_beam_1.keys()[0]
                model_b_2 = amp_beam_2.keys()[0]
                amp_beam_1_chunk = amp_beam_1[model_b_1].data
                amp_beam_2_chunk = amp_beam_2[model_b_2].data

                indata = (
                  sky_model[:,:,wn_chunk],
                  fts_wn_truncated[wn_chunk],
                  spatial_axis,
                  m1_area,
                  td0L,
                  td0R,
                  amp_beam_1_chunk[:,:,wn_chunk],
                  amp_beam_2_chunk[:,:,wn_chunk],
                  beam_angle, 
                  pointing_offset_1,
                  (0, 0),
                  obs_timeline_chunk,
                  True)

                job_id = (times[0], times[-1], beam_angle, pointing_offset_1, wn_chunk.start,
                  wn_chunk.stop)
                job = self.job_server.submit(calculate_visibility,
                  indata, (),
                  ('numpy','collections','matplotlib.pyplot','time','scipy.interpolate',))
                #print '....starting calculation chunk', job_id, len(times)
                jobs.append((job_id, job))

                if len(jobs) >= ncpus:
                    # collect results

                    job_id,job = jobs.popleft()
                    if job is None:
                        raise Exception, 'calculate_visibility has failed for planes %s' % str(job)
                    #print '....calculation chunk', job_id, 'has completed'
                    #print 'unpack', time.time()
                    powers, i1, i2 = job()

                    result_times = powers.keys()
                    for t in result_times:
                        spectra[t].update(powers[t])
                        i1_spectra[t].update(i1[t])
                        i2_spectra[t].update(i2[t])
                    #print 'end unpack', time.time()

            # collect results from remaining jobs
            while len(jobs) > 0:
                job_id,job = jobs.popleft()
                if job is None:
                    raise Exception, 'calculate_visibility has failed for planes %s' % str(job)
                #print '....calculation chunk', job_id, 'has completed'
                # unpacking the result takes a long time - ? tens of seconds
                powers, i1, i2 = job()

                result_times = powers.keys()
                for t in result_times:
                    spectra[t].update(powers[t])
                    i1_spectra[t].update(i1[t])
                    i2_spectra[t].update(i2[t])

            # calculate power for each time
            #print '..assembling spectrum at each timestamp, calculating interferogram value'
            for t in result_times:
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

