Cookbook for PyFIIns
====================
This document describes how to use the FISICA simulator, with some examples.

Process Overview
----------------
The simulation process can be broken into 3 stages.

1. Generate a datacube of the target to be observed. The datacube
is stored in FITS format, the details are described into the companion document
CUBEFORMAT.md. The datacube can be derived from existing observational
data of the object, or built from scratch using the 'cube toolkit' 
described below.

2. Specify the interferometer configuration during the observation 
in the 'instrument' Excel file. This file comprises several sheets, each
describing a particular aspect of the observation. Most users will just want 
to set the waveband of the simulation and the uv pattern to be used.
  
3. Run the simulator, giving the 'instrument' file and the
target datacube as inputs. The simulator will output an html 'weblog' that
can be examined with a browser, and a FITS file containing
the generated interferograms.

4. Run PyDataProcessing with the 'interferograms' file as input. This
will push the data through a data reduction recipe and generate a 'dirty'
cube. In an ideal world this would equal the target cube convolved with
the dirty beam of the interferometer, in reality noise and systematic errors
will further affect the result.

Target Generation - the 'cube toolkit'
--------------------------------------
The simulator can read any target cube whose structure conforms
to CUBEFORMAT.md. The 'makecube' module in the code collection has been written
to make it easier to construct such cubes. Cube construction has four 
parts:  

1. Generate a 2-d image of the target. This can be done by reading
a FITS file with the observed image of a target and interpolating it to the 
required size and spatial resolution (*makecube.MakeImage*). Alternatively, 
the image can be constructed from scratch (*makecube.MakeModelThinRing*,
*MakeModelThickRing* or *MakeModelComet*). The scratch images
are derived from simple models but have the advantage that they are noiseless.

2. Generate the target spectrum (*makecube.MakeSpectrum*). This method 
constructs a greybody spectrum of given temperature, emissivity index and 
peak flux. Some spectral features can be added to the greybody, currently 
available are: 

  * 'forsterite'- a complex feature near 70&mu;m
  * 'protostar' - a set of lines: [OI] 63&mu;m, several water and CO 
lines between 66 and 325&mu;m
  * 'comet' - water lines at 66 and 300&mu;m  

3. Combine the spatial image and the spectrum to produce a cube, write 
the result to a FITS file (*makecube.MakeCube*). The integrated 
flux of the cube equals the input spectrum.

4. Combine cubes. It is possible to build complex targets by adding
together cubes (*makecube.AddCubes*). 
No check is made that the cubes share the same dimensions, if they do not
then a Python error will occur when the add is attempted. The output cube 
inherits its header keywords from the first FITS cube in the 
'in_cubes' list.

The 'Instrument' File
---------------------
The configuration and timeline of the interferometer during the simulated
observation are specified in an Excel spreadsheet. The file contains 
several sheets, each associated with one aspect of the instrument. The file 
'fisica/excel/strawman.xlsx' contains parameters that describe the FISICA
strawman instrument. It should be copied to the working directory and 
modified there.

Most users will only need to edit 2 parts of the strawman file.
First, edit FTSpectrograph/Selection to place a 1 in the Band to be 
simulated, 0 for the others. Second, edit Interferometer/Select in the same
way to specify which uv pattern the instrument is to use in the 
observation. 

Example 1. A simulated observation of protostellar core
--------------------------------------------------------
1. The construction of the
target cube involves the making of 4 separate 'cores', using 
*makecube.MakeModelThickRing* 4 times with different input parameters, then
combining the 4 individual cubes with *makecube.AddCubes*. Thus:

  > In [1]: import makecube  
  >
  > In [2]: i1 = makecube.MakeModelThickRing(rinner=-0.4, router=1.0,  
  >	tilt=20.0, rot=-90.0, xcentre=-0.7, ycentre=0.2, max_baseline=80.0,  
  >	wn_min=100.0, wn_max=200.0)  
  >
  > In [3]: image1 = i1.run()
  >
  > In [4]: s1 = makecube.MakeSpectrum(temperature=37, beta=-1.3, 
  >   wn_min=100.0, wn_max=200.0, wn_step=0.5, peak_flux=10.0, 
  >   spectral_features=['forsterite'])  
  >  
  > In [5]: spectrum1 = s1.run()  
  >
  > In [6]: c1 = makecube.MakeCube(image=image1, spectrum=spectrum1,  
  >   cubename='cube1.fits')  
  >
  > In [7]: c1.run()  
  >
  >... repeat this process n times to produce n 'cores' in n cubes, each  
  >... with different parameters. This will give you n FITS files: 'cube1.fits'  
  >... to 'cuben.fits'.  
  >
  > In [20]: makecube.AddCubes(in_cubes=['cube1.fits', 'cube2.fits',  
  >..., 'cuben.fits'])

2. Specify the observation details by copying 'strawman.xlsx' to your 
working directory and modifying it as follows:
  * In FTSpectrograph select Band 2 (100-200cm-1).
  * In Interferometer select 'Const L spiral' as the uv pattern.

3. Run the simulator:
  > In [1]: import pyfiins  
  >
  > In [2]: f = pyfiins.PyFIInS(instrument_spreadsheet='strawman_copy.xlsx',
  >   sky_file='your_target.fits',
  >   beam_model_dir='/Users/jfl/test2beam/fisica-pbmodels/beam_models/Smooth_Walled_Horn/Band 4 GRASP')  

This will produce a weblog html tree rooted at 'fisica-sim-yyyymmddThhmmss'
dated at the time of the run, and a matching FITS file with the 
simulated interferograms 'yyyymmddThhmmss.fits'

3. Reduce the interferograms.
  > In [3]: import pydataprocessing  
  >
  > In [4]: d=pydataprocessing.PyDataProcessing('your_interferograms.fits',  
  >   data_quality='pure')  
  >
  > In [5]: d.reduce()

The 'data_quality' parameter in the PyDataProcessing constructor  
specifies whether the 'pure' or 'noisy' interferograms are to be used.
The pure interferograms are the direct output from the interferometric
simulation, the noisy ones have detector noise, background noise, 
cosmic rays and 1/f noise added.


Example 2. A simulated observation of the accretion disk around a young star
----------------------------------------------------------------------------
TBD 



