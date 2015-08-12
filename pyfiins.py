from __future__ import absolute_import

import astropy.io.fits as pyfits
import collections
import datetime
import pp

import addnoise
import detectornoise
import fts
import loadparameters
import noise
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
    def __init__(self, sky_file='SkySparams.xlsx', sky_sheet='1point',
      instrument_spreadsheet='FIInS_Instrument_cor3.xlsx',
      beam_model_dir='.'):
        # is this an Excel spreadsheet with source parameters to simulate?
        # or is it a FITS file containing a cube simulated beforehand?
        self.sky_spreadsheet = None
        self.sky_sheet = None
        self.sky_fits = None
        try:
            book = xlrd.open_workbook(sky_file)
            self.sky_spreadsheet = sky_file
            self.sky_sheet = sky_sheet
            print 'sky simulation parameters will be read from Excel file %s' % \
              sky_file
        except:
            pass
        try:
            hdulist = pyfits.open(sky_file)
            self.sky_fits = sky_file
            print 'a simulated sky datacube will be read from FITS file %s' % \
              sky_file
        except:
            pass

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
          sky_spreadsheet=self.sky_spreadsheet,
          sky_sheet=self.sky_sheet,
          sky_fits=self.sky_fits,
          instrument_spreadsheet=self.instrument_spreadsheet)
        obs_specification = loadparams.run()
        self.result['loadparameters'] = obs_specification

        # generate information on the FTS 
        ftsd = fts.FTS(parameters=obs_specification)
        self.result['fts'] = ftsd.run()
        print ftsd

        # generate information on the flux collectors 
        tel = telescope.Telescope(parameters=obs_specification)
        self.result['telescope'] = tel.run()
        print tel
       
        # generate UV map
        uvmapgen = uvmapgenerator.UVMapGenerator(
          parameters=obs_specification,
          previous_results=self.result)
        self.result['uvmapgenerator'] = uvmapgen.run()
        print uvmapgen

        # calculate noise
        rumble = noise.Noise(parameters=obs_specification,
          previous_results=self.result)
        self.result['noise'] = rumble.run()
        print rumble

        # calculate detector noise
        dn = detectornoise.DetectorNoise(parameters=obs_specification,
          previous_results=self.result)
        self.result['detectornoise'] = dn.run()
        print dn

        # construct sky
        if self.sky_spreadsheet is not None:
            skygen = skygenerator.SkyGenerator(
              parameters=obs_specification,
              previous_results=self.result)
            self.result['skymodel'] = skygen.run()
            print skygen
        elif self.sky_fits is not None:
            skyload = skyloader.SkyLoader(
              parameters=obs_specification,
              previous_results=self.result)
            self.result['skymodel'] = skyload.run()
            print skyload
   
        # generate primary beams
        primarybeamsgen = pbmodelgenerator.PrimaryBeamsGenerator(
          previous_results=self.result,
          beam_model_dir = self.beam_model_dir,
          job_server=self.job_server)
        self.result['primarybeams'] = primarybeamsgen.run()
        print primarybeamsgen
   
        # generate observation framework
        timeline = timelinegenerator.TimeLineGenerator(
          previous_results=self.result)
        self.result['timeline'] = timeline.run()
        print timeline
   
        # calculate interferograms
        obs = observe.Observe(
          parameters=obs_specification,
          previous_results=self.result,
          job_server=self.job_server)
        self.result['observe'] = obs.run()
        print obs

        # add noise, cosmic rays, detector time constant
        with_errors = addnoise.AddNoise(
          parameters=obs_specification,
          previous_results=self.result)
        self.result['addnoise'] = with_errors.run()
        print with_errors

        # write out the interferograms as FITS files
        fits = writefits.WriteFITS(previous_results=self.result)
        self.result['writefits'] = fits.run()   
        print fits

        # construct html description of result
        htmlrenderer = renderer.Renderer(result=self.result)
        htmlrenderer.run(prefix='sim')

    def __repr__(self):
        return 'FISICA simulator'

