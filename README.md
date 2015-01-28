# fisica
pyfiins - a Python tool for simulating a Double Fourier Interferometer in space. 

The FISICA Project
------------------
FISICA stands for 'Far Infra-red Space Interferometer Critical Assessment': 
Scientific Definition and Technology Development for the Next Generation 
THz Space Interferometer

It is funded by the European Union as part of the Seventh Framework Programme 
- FP7. The project is funded for three years, starting January 2013.

pyfiins
-------
This is a python program being written to simulate the data taken by
a 'strawman design' double Fourier interferometer in space.

Currently, it is run from IPython (see ipython.org) built on python 2.7 on a 
Mac.

Installation
------------
This describes the installation procedure for a Mac. We will attempt to automate
the installation process in the near future and make sure that it works for
other OS.

1. Follow ipython.org/install instructions to install Anaconda and then use 
that to install IPython and the other modules needed:

 numpy (for fast array processing)
 scipy (scientific functions and constants)
 matplotlib (for plotting) 
 parallelpython (for parallel processing)
 
Yo will also need to install:
 mako and twitter bootstrap (weblog generation)
 
If the installation has worked then you should be able to type ipython, then:

In [1]: import numpy
In [2]: import scipy
In [3]: import matplotlib
In [4]: import pp
In [5]: import mako

with no errors.

2. Check out the trunk code from GitHub.

3. Copy the .xlsx files from fisica/excel into the directory you 
will be running python from - that is where the code will look for them.

4. Start up python then:

In [1]: import sys
  tell python where to look for the code
In [2]: sys.path.append('/users/jfl/fisica')

In [3]: import fiins
  construct the simulator specifying a file describing the sky model 
  to use
In [4]: f=fiins.Fiins(sky_spreadsheet='SkyTest.xlsx', sky_sheet='Master')
In [5]: f.simulate()

It should run for several minutes, producing lots of debug statements 
and a few 'warnings', eventually finishing with:

.....
rendering observe
rendering reduceinterferogram
rendering dirtyimage
rendering cleanimage

If it has worked then you should see a new directory in your working 
directory with name of form 'fisisca-yyyymmddThhmmss', which contains 
the browsable results for the run just made. Point your browser to 
that directory, scroll down to the bottom and click on 
'uvmapgenerator.html' whereupon the navigable weblog display should appear.
