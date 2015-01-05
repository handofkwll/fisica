from __future__ import absolute_import

import collections
import pp

import cleanimage
import cubeparameters
import dirtyimage
import fts
import loadparameters
import observe
import primarybeamsgenerator
import renderer
import skygenerator
import synthesisbeamsgenerator
import timelinegenerator
import uvmapgenerator
import uvspectra


class pyfiins(object):
    """FISICA simulator.
    """
    def __init__(self, sky_spreadsheet='SkySparams.xlsx', sky_sheet='1point'):
        self.sky_spreadsheet = sky_spreadsheet
        self.sky_sheet = sky_sheet
        self.result = collections.OrderedDict()

        # start parallel python (pp), find the number of CPUS available
        ppservers = ()
        self.job_server = pp.Server(ppservers=ppservers)
        print 'pyfiins starting pp with %s workers' % \
          self.job_server.get_ncpus()

    def simulate(self):
        # read parameters
        loadparams = loadparameters.LoadParameters(
          sky_spreadsheet=self.sky_spreadsheet,
          sky_sheet=self.sky_sheet)
        obs_specification = loadparams.run()
        self.result['loadparameters'] = obs_specification

        # generate information on the FTS 
        ftsd = fts.FTS(parameters=obs_specification)
        self.result['fts'] = ftsd.run()
        print ftsd

        # generate UV map
        uvmapgen = uvmapgenerator.UVMapGenerator(parameters=obs_specification,
          previous_results=self.result)
        self.result['uvmapgenerator'] = uvmapgen.run()
        print uvmapgen

        # generate parameters of cube to be simulated
        cubeparams = cubeparameters.CubeParameters(
          previous_results=self.result)
        self.result['cubeparameters'] = cubeparams.run()
        print cubeparams

        # generate synthesis beams - dirty and clean
        synthbeamsgen = synthesisbeamsgenerator.SynthesisBeamsGenerator(
          previous_results=self.result, job_server=self.job_server)
        self.result['synthesisbeams'] = synthbeamsgen.run()
        print synthbeamsgen

        # generate primary beams
        primarybeamsgen = primarybeamsgenerator.PrimaryBeamsGenerator(
          previous_results=self.result, job_server=self.job_server)
        self.result['primarybeams'] = primarybeamsgen.run()
        print primarybeamsgen

        # construct sky
        skygen = skygenerator.SkyGenerator(parameters=obs_specification,
          previous_results=self.result)
        self.result['skygenerator'] = skygen.run()
        print skygen

        # generate observation framework
        timeline = timelinegenerator.TimeLineGenerator(
          previous_results=self.result)
        self.result['timeline'] = timeline.run()
        print timeline

        # calculate interferograms
        obs = observe.Observe(parameters=obs_specification,
          previous_results=self.result, job_server=self.job_server)
        self.result['observe'] = obs.run()  
        print obs

        # recover spectra from interferograms
#        uvspec = uvspectra.UVspectra(previous_results=self.result,
#          job_server=self.job_server)
#        self.result['uvspectra'] = uvspec.run()  
#        print uvspec

        # construct dirty image
#        dirty = dirtyimage.DirtyImage(previous_results=self.result,
#          job_server=self.job_server)
#        self.result['dirtyimage'] = dirty.run()  
#        print dirty

#        # construct clean image
#        clean = cleanimage.CleanImage(previous_results=self.result,
#          job_server=self.job_server)
#        self.result['cleanimage'] = clean.run()  
#        print clean

        # construct html description of result
        self.render()
        
    def import_result(self):
        print 'does nothing'

    def render(self):
        htmlrenderer = renderer.Renderer(result=self.result)
        htmlrenderer.run()

    def __repr__(self):
        return 'FISICA simulator'

