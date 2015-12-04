from __future__ import absolute_import

import astropy.io.fits as pyfits
import collections
import datetime
import pp

import addnoise
import backgroundnoise
import detectornoise
import fts
import loadparameters
import observe_interpolate as observe
import pbmodelgenerator
import renderer
import skygenerator
import skyloader
import telescope
import timelinegenerator
import uvmapgenerator
import writefits
import xlrd


class PyFIInS(object):
    """FISICA simulator.
    """
    def __init__(self, sky_file='ringcube.fits',
      instrument_spreadsheet='FIInS_Instrument.xlsx',
      beam_model_dir='.'):
        self.sky_fits = sky_file
        self.instrument_spreadsheet = instrument_spreadsheet
        self.beam_model_dir = beam_model_dir

        self.result = collections.OrderedDict()

        # start parallel python (pp), find the number of CPUS available
        ppservers = ()
        self.job_server = pp.Server(ppservers=ppservers)
        print 'PyFIInS starting pp with %s workers' % \
          self.job_server.get_ncpus()

    def simulate(self):
        # store data and time of run
        now = datetime.datetime.today()
        self.result['runtime'] = now.strftime('%Y%m%dT%H%M%S')

        # read parameters
        loadparams = loadparameters.LoadParameters(
          sky_fits=self.sky_fits,
          instrument_spreadsheet=self.instrument_spreadsheet)
        obs_specification = loadparams.run()
        self.result['loadparameters'] = obs_specification
        del loadparams

        # generate information on the FTS 
        ftsd = fts.FTS(parameters=obs_specification)
        self.result['fts'] = ftsd.run()
        print ftsd

        # generate information on the flux collectors 
        tel = telescope.Telescope(parameters=obs_specification)
        self.result['telescope'] = tel.run()
        print tel
        del tel
       
        # generate UV map
        uvmapgen = uvmapgenerator.UVMapGenerator(
          parameters=obs_specification,
          previous_results=self.result)
        self.result['uvmapgenerator'] = uvmapgen.run()
        print uvmapgen
        del uvmapgen

        # calculate background noise
        background = backgroundnoise.BackgroundNoise(
          parameters=obs_specification, previous_results=self.result)
        self.result['backgroundnoise'] = background.run()
        print background
        del background

        # construct sky
        skyload = skyloader.SkyLoader(
          parameters=obs_specification,
          previous_results=self.result)
        self.result['skymodel'] = skyload.run()
        print skyload
        del skyload   

        # generate primary beams
        primarybeamsgen = pbmodelgenerator.PrimaryBeamsGenerator(
          previous_results=self.result,
          beam_model_dir = self.beam_model_dir,
          job_server=self.job_server)
        self.result['primarybeams'] = primarybeamsgen.run()
        print primarybeamsgen
        del primarybeamsgen   

        # generate observation framework
        timeline = timelinegenerator.TimeLineGenerator(
          previous_results=self.result)
        self.result['timeline'] = timeline.run()
        print timeline
        del timeline

        # calculate detector noise
        dn = detectornoise.DetectorNoise(parameters=obs_specification,
          previous_results=self.result)
        self.result['detectornoise'] = dn.run()
        print dn
        del dn   

        # calculate interferograms
        obs = observe.Observe(
          parameters=obs_specification,
          previous_results=self.result,
          job_server=self.job_server)
        self.result['observe'] = obs.run()
        print obs
        del obs

        # add noise, cosmic rays, detector time constant
        with_errors = addnoise.AddNoise(
          parameters=obs_specification,
          previous_results=self.result)
        self.result['addnoise'] = with_errors.run()
        print with_errors
        del with_errors

        # write out the interferograms as FITS files
        fits = writefits.WriteFITS(previous_results=self.result)
        self.result['writefits'] = fits.run()   
        print fits
        del fits

        # construct html description of result
        htmlrenderer = renderer.Renderer(result=self.result)
        htmlrenderer.run(prefix='sim')

    def __repr__(self):
        return 'FISICA simulator'

