"""This module contains classes and methods used to save data from the
simulated observation to a FITS-IDI file. 
"""

from __future__ import absolute_import

import datetime
import ephem
import numpy as np
import os.path
import astropy.constants as constants
import astropy.io.fits as pyfits
import astropy.time as astro_time

# FITS IDI python module imports
from pyFitsidi import *
from astroCoords import *

def config_antenna(tbl):
  """ Configures the antenna table.
  
  Parameters
  ----------
  tbl: pyfits.hdu
    table to be configured
  """

  antenna = tbl.data
  for i in range(0,tbl.data.size):
    antenna[i]['ANNAME']      = 'FIRI_%i' % i
    antenna[i]['ANTENNA_NO']  = i
    antenna[i]['ARRAY']       = 1
    antenna[i]['FREQID']      = 1
    antenna[i]['NO_LEVELS']   = 12
    antenna[i]['POLTYA']      = 'X'
    antenna[i]['POLTYB']      = 'X'
    antenna[i]['POLAA']       = np.array(0)
    antenna[i]['POLAB']       = np.array(0)

  tbl.data = antenna
  return tbl

def config_array_geometry(tbl, antenna_array):
  """  Configures the array_geometry table with FIRI values

  Parameters
  ----------
  tbl: pyfits.hdu
    table to be configured
  antenna_array: numpy.array
    an array of xyz coordinates of the antenna locations (offsets) in METERS
    from the array centre (this is a keyword in the header unit)
    e.g. 
  """
  
  geometry = tbl.data

  # X-Y-Z in metres
  xyz_m = antenna_array

  for i in range(0,tbl.data.size):
    geometry[i]['ANNAME']  = 'COLL_%i' % i
    geometry[i]['STABXYZ'] = xyz_m[i]
    geometry[i]['DERXYZ']  =  0
    #geometry[i]['ORBPARM'] = 0
    # The antenna 'number' NOSTA is 1-based
    geometry[i]['NOSTA']   = i + 1
    geometry[i]['MNTSTA']  = 2
    # NOTE: Aperture arrays are given code 6, but not supported by CASA
    geometry[i]['STAXOF']  = np.array([0,0,0])
    geometry[i]['DIAMETER'] = 2.0

  tbl.data = geometry

  return tbl

def config_frequency(tbl, external_params):
  """
  Configures the frequency table.
  
  Parameters
  ----------
  tbl: pyfits.hdu
    table to be configured
  """

  frequency = tbl.data[0]

  frequency['FREQID']         = 1
  frequency['BANDFREQ']       = 0         # This is offset from REF_FREQ, so zero!
  frequency['CH_WIDTH']       = external_params['CHAN_BW']
  frequency['TOTAL_BANDWIDTH']= external_params['NCHAN'] * \
    external_params['CHAN_BW']
  frequency['SIDEBAND']       = 1

  tbl.data[0] = frequency

  return tbl

def config_source(tbl, source):
  """  Configures the source table.
  
  Parameters
  ----------
  tbl: pyfits.hdu
    table to be configured
  source: ephem.fixedBody
    source to be phased to (use makeSource())
  """
  
  # Stupidly using source as a variable name twice
  source_ra   = np.rad2deg(source._ra)
  source_dec  = np.rad2deg(source._dec)
  source_name = source.name
  
  print('Source is: %s'%source.name)
  
  source = tbl.data[0]
  
  source['SOURCE_ID'] = 1
  source['SOURCE']    = source_name
  source['VELDEF']    = 'RADIO'
  source['VELTYP']    = 'GEOCENTR'
  source['FREQID']    = 1
  source['RAEPO']     = source_ra
  source['DECEPO']    = source_dec
  source['EQUINOX']   = 'J2000'
  
  # Things I'm just making up
  source['IFLUX']    = 0
  source['QFLUX']    = 0
  source['UFLUX']    = 0
  source['VFLUX']    = 0
  source['ALPHA']    = 0
  source['FREQOFF']  = 0
  
  tbl.data[0] = source
  
  return tbl


class Array(ephem.Observer):
    """ An antenna array class.
    
        Based on pyEphem's Observer class.
        Probably very similar to the one in AIPY.
        
        Parameters
        ----------
        lat: dd:mm:ss
          latitude of array centre, e.g. 44:31:24.88
        long: dd:mm:ss
          longitude of array centre, e.g. 11:38:45.56
        elev: float
          elevation in metres of array centrem e.g. 28.0
        date: datetime object
          date and time of observation, e.g. datetime.now()
        antennas: np.array([x,y,z])
          numpy array of antenna positions, in xyz coordinates in meters,
          relative to the array centre.
        
    """

    def __init__(self, lat, long, elev, date, antennas):
        super(Array, self).__init__()
        self.lat = lat
        self.long = long
        self.elev = elev
        self.date = date
        self.antennas = antennas

    def update(self, date):
        """Update antenna with a new datetime"""
        self.date = date


class WriteFITS_IDI(object):
    """Class to write out observed uv data to FITS-IDI File.

    Contains methods:
    __init__
    run
    __repr__
    """

    def __init__(self, previous_results):
        """Constructor.

        Keyword parameters:
        previous_results - simulator result structure.
        """
        self.previous_results = previous_results
        self.result = {}

    def run(self):
        """Method to write the FITS data.
        """
        #print 'WriteFITS_IDI.run'

        # construct the name of the file
        readfits = self.previous_results['readfits']
        obs_date = readfits['obs date']
        idifitsfile = '%s.idi.fits' % obs_date

        configxml = 'firi.xml'

        # midnight on date to Julian day
        obs_date_midnight = astro_time.Time('%s-%s-%sT00:00:00' %
          (obs_date[:4], obs_date[4:6], obs_date[6:8]), format='isot')
        obs_date_midnight = obs_date_midnight.jd

        rdate = astro_time.Time(obs_date_midnight, format='jd',
          out_subfmt='date')
        rdate = rdate.iso

        # number of days after midnight at obs start
        obs_date_time = astro_time.Time('%s-%s-%s:%s:%s' %
          (obs_date[:4], obs_date[4:6], obs_date[6:11], obs_date[11:13],
          obs_date[13:]), format='isot')
        obs_date_time = obs_date_time.jd - obs_date_midnight

        # get specific items from the results that will be need in
        # the reduction
        reduce_interferogram = self.previous_results['reduceinterferogram']
        data_quality = reduce_interferogram['data_quality']
        scan_uvspectra = reduce_interferogram['scan_uvspectra']

        wavenumber = scan_uvspectra[0].wavenumber

        # construct lists of the values to be stored in each Table column
        n_uvspectra = max(scan_uvspectra.keys()) + 1
        mcomplex = 3
        mstokes = 1
        mfreq = len(wavenumber)
        mra = 1
        mdec = 1

        uv_data = np.zeros([n_uvspectra, mdec, mra, mfreq, mstokes, mcomplex])
        u = np.zeros([n_uvspectra])
        v = np.zeros([n_uvspectra])
        w = np.zeros([n_uvspectra])
        dates = np.zeros([n_uvspectra])
        times = np.zeros([n_uvspectra])
        baselines = np.zeros([n_uvspectra], dtype=np.int)
        freqid = np.ones([n_uvspectra], dtype=np.int)

        for k,val in scan_uvspectra.items():
            uv_data[k,0,0,:,0,0] = val.spectrum.real
            uv_data[k,0,0,:,0,1] = val.spectrum.imag
            uv_data[k,0,0,:,0,2] = np.ones(val.spectrum.real.shape)
            u[k] = np.mean(val.baseline_x)
            v[k] = np.mean(val.baseline_y)
            w[k] = np.mean(val.baseline_z)
            dates[k] = obs_date_midnight
            times[k] = obs_date_time + (np.mean(val.time) / (3600 * 24))
            baselines[k] = 258

        # external_params is referred to inside config.xml and can be
        # used to set parameters there
        light_speed = constants.c.to('m/s').value
        external_params = {'NCHAN':len(wavenumber),
                           'RDATE':rdate,
                           'REF_FREQ':0.0 * 100 * light_speed,
                           'CHAN_BW':np.abs(wavenumber[1] - wavenumber[0]) * \
                             100 * light_speed}

        print "Out: %s\nConfig: %s"%(idifitsfile, configxml)

        print('\nConfiguring Array geography')
        print('--------------------------')
        # Meaningless numbers, hopefully not needed by any CASA method 
        # that we want to use
        (latitude, longitude, elevation) = ('00:00:00.00', '00:00:00.00', 0)
        now = datetime.datetime.now()

        # Make ourselves an Array (pyEphem observer)
        array_geometry_m = np.array([
          [0.0, 0.0, 0.0],
          [0.0, 80.0, 0.0]], dtype = 'float32')
        beach = Array(lat=latitude, long=longitude, elev=elevation, date=now,
          antennas=array_geometry_m)

        print('\nConfiguring phase source')
        print('--------------------------')
        # The source is our phase centre for UVW coordinates
        line = "%s,f,%s,%s,%s,%d" % ('Deep Space', '00:00:00',
          '00:00:00', '1', 2000)
        source = ephem.readdb(line)
        source.compute(beach)
        print "Name: %s \nRA: %s \nDEC: %s"%(source.name, source.ra, source.dec)

        # Make a new blank FITS HDU
        print('\nCreating PRIMARY HDU')
        print('------------------------------------')
        hdu = make_primary(config=configxml, external_params=external_params)
        print repr(hdu.header)

        # Go through and generate required tables
        print('\nCreating ARRAY_GEOMETRY')
        print('------------------------------------')
        tbl_array_geometry = make_array_geometry(config=configxml, num_rows=2,
          external_params=external_params)
        tbl_array_geometry = config_array_geometry(tbl_array_geometry,
          array_geometry_m)
        print repr(tbl_array_geometry.header)

        print('\nCreating FREQUENCY')
        print('------------------------------------')
        tbl_frequency = make_frequency(config=configxml, num_rows=1,
          external_params=external_params)
        tbl_frequency = config_frequency(tbl_frequency,
          external_params=external_params)
        print repr(tbl_frequency.header)

        print('\nCreating SOURCE')
        print('------------------------------------')
        tbl_source = make_source(config=configxml, num_rows=1,
          external_params=external_params)
        tbl_source = config_source(tbl_source, source)
        print repr(tbl_source.header)

        print('\nCreating ANTENNA')
        print('------------------------------------')
        tbl_antenna = make_antenna(config=configxml, num_rows=2,
          external_params=external_params)
        tbl_antenna = config_antenna(tbl_antenna)
        print repr(tbl_antenna.header)

        print('\nCreating UV_DATA')
        print('------------------------------------')

        print 'Data dimensions: %i dumps, %i chans, %i pols, %i data' % (
          n_uvspectra, mfreq, mstokes, mcomplex)

        print('Generating blank UV_DATA rows...')
        tbl_uv_data = make_uv_data(config=configxml, num_rows=n_uvspectra,
          external_params=external_params)

        timesorted = np.argsort(times)

        for k in timesorted:
            tbl_uv_data.data[k]['FLUX'] = uv_data[k,0,0,:,0,:].ravel()
            tbl_uv_data.data[k]['UU'] = u[k] / light_speed
            tbl_uv_data.data[k]['VV'] = v[k] / light_speed
            tbl_uv_data.data[k]['WW'] = w[k] / light_speed
            tbl_uv_data.data[k]['BASELINE'] = baselines[k]
            tbl_uv_data.data[k]['DATE'] = dates[k]
            tbl_uv_data.data[k]['TIME'] = times[k]
            tbl_uv_data.data[k]['SOURCE'] = 1
            tbl_uv_data.data[k]['FREQID'] = 1
            tbl_uv_data.data[k]['INTTIM'] = 3

        print repr(tbl_uv_data.header)
  
        hdulist = pyfits.HDUList(hdus=
           [hdu,
           tbl_array_geometry,
           tbl_source, 
           tbl_frequency,
           tbl_antenna,
           tbl_uv_data])

        print('Verifying integrity...')            
        hdulist.verify()
  
        if(os.path.isfile(idifitsfile)):
            print('Removing existing file...')
            os.remove(idifitsfile)
        print('Writing to file...')
        hdulist.writeto(idifitsfile)

        print('Done.')

        self.result['idifitsfile'] = idifitsfile

        return self.result

    def __repr__(self):
        return '''
WriteFITS_IDI:
  FITS file : {fitsfile}
'''.format(fitsfile=self.result['idifitsfile'])
