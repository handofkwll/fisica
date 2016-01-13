"""This module contains the classes and methods used to calculate the
background noise falling on the FIRI detector.
"""

from __future__ import absolute_import
__author__ = 'rjuanola'
import collections
import numpy as np
import scipy.constants as sc

def black_body(temperature, wavenumber):
    """Function to calculate Planck function.

    Keyword arguments:
    temperature - Kelvin
    wavenumber  - frequency [cm-1]

    Returns:
    value of Planck function
    """
    freq = wavenumber * sc.c * 100
    bb = 2 * sc.h * freq**3 / \
      (sc.c**2 * (np.exp((sc.h * freq) / (sc.k * temperature)) - 1.0))
    return bb

def NEP2(eta_feed, throughput, emissivity, trans2det, temp, freq):
    """Function to calculate Noise Equivalent Power (NEP) following the
    method described in Roser Juanola's thesis p117, which references
    J.M.Lamarre, "Photon noise in photometric instruments at far-infrared
    and sub-millimeter wavelengths", Appl. Opt., 25, 870-876,
    Mar 1986.

    Keyword arguments:
    eta_feed   -- detector efficiency
    throughput -- A * omega
    emissivity -- emissivity of emitter
    trans2det  -- transmission to detector
    temp       -- emitter temperature [K]
    freq       -- frequency grid [Hz]

    Returns:
    value of NEP**2
    """
    den = np.exp(sc.h * freq / (sc.k * temp)) - 1
    a = 1 + (emissivity * trans2det * eta_feed / den)
    b = throughput * emissivity * trans2det * eta_feed * freq**4 / den
    cte = 4 * (sc.h / sc.c)**2
    d = a * b * (freq[2] - freq[1])
    nep2 = cte * np.sum(d)
    return nep2


class BackgroundNoise(object):
    """Class to generate the NEP of the simulated observation.

    Contains methods:
    __init__
    run
    __repr__
    """

    def __init__(self, parameters, previous_results):
        """Constructor.
        Parameters:
          parameters       - Structure containing parameters defining
                             the simulation, read from Excel.
          previous_results - Current results structure of the simulation run.
        """
        self.parameters = parameters
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        """Method invoked to do the work.
        """

        # Loading background parameters
        background = self.parameters['substages']['Background']
        row = background['CMB temp'].keys()[0]
        Tcmb = background['CMB temp'][row]
        row = background['CIB temp'].keys()[0]
        Tcib = background['CIB temp'][row]
        row = background['CIB emissivity @ 200um'].keys()[0]
        Ecib = background['CIB emissivity @ 200um'][row]
        row = background['Zodi dust temp'].keys()[0]
        Tzodi = background['Zodi dust temp'][row]
        row = background['Zodi emissivity @ 30um'].keys()[0]
        Ezodi = background['Zodi emissivity @ 30um'][row]
        row = background['Solar temp'].keys()[0]
        Tsun = background['Solar temp'][row]
        row = background['Zodi scattering coeff'].keys()[0]
        Szodi = background['Zodi scattering coeff'][row]
        row = background['Tstray'].keys()[0]
        Tstray = background['Tstray'][row]
        row = background['Stray coeff'].keys()[0]
        stray_coeff = background['Stray coeff'][row]

        #Loading cold optics parameters
        coldoptics = self.parameters['substages']['ColdOptics']
        row = coldoptics['Optical efficiency (Oeff)'].keys()[0]
        Oeff = coldoptics['Optical efficiency (Oeff)'][row]
        row = coldoptics['Number of mirrors (NmirrorsL)'].keys()[0]
        NmirrorsL = coldoptics['Number of mirrors (NmirrorsL)'][row]
        row = coldoptics['Cold optics box temp (Tbox)'].keys()[0]
        Tbox = coldoptics['Cold optics box temp (Tbox)'][row]
        row = coldoptics['Instrument emissivity (InstEmi)'].keys()[0]
        InstEmi = coldoptics['Instrument emissivity (InstEmi)'][row]
        row = coldoptics['Cold optics mirror emissivity (COEmis)'].keys()[0]
        COEmis = coldoptics['Cold optics mirror emissivity (COEmis)'][row]
        row = coldoptics['Dichroic transmission (DicT)'].keys()[0]
        DicT = coldoptics['Dichroic transmission (DicT)'][row]
        row = coldoptics['Combiner T (CoT)'].keys()[0]
        CoT = coldoptics['Combiner T (CoT)'][row]
        row = coldoptics['Dichroic reflectivity (DicR)'].keys()[0]
        DicR = coldoptics['Dichroic reflectivity (DicR)'][row]
        row = coldoptics['Combiner T (CoT)'].keys()[0]
        CoT = coldoptics['Combiner T (CoT)'][row]
        row = coldoptics['Combiner R (CoR)'].keys()[0]
        CoR = coldoptics['Combiner R (CoR)'][row]
        row = coldoptics['Window efficiency (Weff)'].keys()[0]
        Weff = coldoptics['Window efficiency (Weff)'][row]
        row = coldoptics['NIRFIR dichroic Tx (NFdicT)'].keys()[0]
        NFdicT = coldoptics['NIRFIR dichroic Tx (NFdicT)'][row]
        row = coldoptics['Filter Tx (FilT)'].keys()[0]
        FilT = coldoptics['Filter Tx (FilT)'][row]

        # L and R refer to the lights paths from the left and right
        # collectors. In the FIRI design described in Roser's thesis
        # one arm needs an extra mirror to direct the beam onto the 
        # beamsplitter
        NmirrorsR = NmirrorsL - 1

        # Loading warm optics parameters
        warmoptics = self.parameters['substages']['WarmOptics']
        row = warmoptics['BSM temperature'].keys()[0]
        Tbsm = warmoptics['BSM temperature'][row]
        row = warmoptics['BSM emissivity'].keys()[0]
        Ebsm = warmoptics['BSM emissivity'][row]
        row = warmoptics['K-mirror temperature'].keys()[0]
        Tk = warmoptics['K-mirror temperature'][row]
        row = warmoptics['K-mirror emissivity'].keys()[0]
        Ek = warmoptics['K-mirror emissivity'][row]
        row = warmoptics['Number of mirrors in k-mirror'].keys()[0]
        Nk = warmoptics['Number of mirrors in k-mirror'][row]

        # get range of wavenumbers being simulated
        fts = self.previous_results['fts']
        fts_wn_truncated = fts['fts_wn_truncated']

        # load collector parameters
        tel = self.previous_results['telescope']
        Ttel = tel['mirror_temp']
        Etel = tel['telescope_emissivity']
        Dtel = tel['m1_diameter']

        # load detector parameters
        detectors = self.parameters['substages']['Detectors']
        row = detectors['Detector optical absorption efficiency'].keys()[0]
        det_eff = detectors['Detector optical absorption efficiency'][row]

        # Cold optics transmission
        # ..Box left arm to detector 1
        # .. windw * NIRdichroic * coldopt**nopt * dichroic * filt * combiner
        t4L = Weff * NFdicT * (1 - COEmis)**NmirrorsL * DicT * FilT * CoT
        # ..Box right arm to detector 1
        t4R = Weff * NFdicT * (1 - COEmis)**NmirrorsR * DicT * FilT * CoR 

        # Cold optics emissivities
        e4L = 1 - t4L
        e4R = 1 - t4R

        # Warm optics transmissions and emissivities
        t1 = 1 - Etel            #telescope primary
        t2 = 1 - Etel            #telescope secondary
        tbsm = 1 - Ebsm          #beam steering mirror
        tk = 1                   #no k mirror (1-Ek)^Nk;          #k mirror, right side
        t3 = 1                   #stray light
        e1 = 1 - t1
        e2 = 1 - t2
        ebsm = 1 - tbsm
        ek = 1 - tk
        e3 = Etel

        #transmissions to detector 1 from left
        td0L = t1 * t2 * tbsm * t4L
        td1L = t2 * tbsm * t4L

        td12L = tbsm * t4L

        td2L = t4L
        td3L = t4L
        td4L = 1

        # transmissions to detector 1 from right
        td0R = t1 * t2 * tbsm * tk * t4R
        td1R = t2 * tbsm * tk * t4R

        td12R = tbsm * tk * t4R
        td13R = tk * t4R

        td2R = t4R
        td3R = t4R
        td4R = 1

        # Computing neps for instrument
        rm = fts_wn_truncated * 100

        # Assuming single mode optics A * omega = etendue = lambda**2
        etendue = (1 / rm)**2

        # NEP due to emission from primary
        NEPL1 = NEP2(det_eff, etendue, e1, td1L, Ttel, rm * sc.c)
        # NEP due to emission from secondary
        NEPL2 = NEP2(det_eff, etendue, e2, td2L, Ttel, rm * sc.c)
        # NEP due to emission from beamsplitter
        NEPL3 = NEP2(det_eff, etendue, ebsm, td12L, Ttel, rm * sc.c)
        # NEP due to stray light not from collector optics, this is a fudge
        NEPL4 = NEP2(det_eff, etendue, e3, stray_coeff, Tstray, rm * sc.c)
        # NEP due to emission from cold box
        NEPL5 = NEP2(det_eff, etendue, e4L, td4L, Tbox, rm * sc.c)

        NEPR1 = NEP2(det_eff, etendue, e1, td1R, Ttel, rm * sc.c)
        NEPR2 = NEP2(det_eff, etendue, e2, td2R, Ttel, rm * sc.c)
        NEPR3 = NEP2(det_eff, etendue, ebsm, td12R, Tbsm, rm * sc.c)
        NEPR4 = NEP2(det_eff, etendue, e3, stray_coeff, Tstray, rm * sc.c)
        NEPR5 = NEP2(det_eff, etendue, e4R, td4R, Tbox, rm * sc.c)

        NEPL = NEPL1 + NEPL2 + NEPL3 + NEPL4 + NEPL5
        NEPR = NEPR1 + NEPR2 + NEPR3 + NEPR4 + NEPR5

        # TOTAL BACKGROUND NEP
        # ..from instrument)
        NEPtotInst2 = (np.sum(NEPL) + np.sum(NEPR))
        NEPtotInst = np.sqrt(NEPtotInst2)

        # ..from CMB
        EcmbL = 1 - e1 - e2 - ebsm - e4L
        EcmbR = 1 - e1 - e2 - ebsm - e4R

        NEPcmbL = NEP2(det_eff, etendue, EcmbL, td3L, Tcmb, rm*sc.c)
        NEPcmbR = NEP2(det_eff, etendue, EcmbR, td3R, Tcmb, rm*sc.c)
        NEPcmb2 = NEPcmbL + NEPcmbR
        NEPcmb = np.sqrt(NEPcmb2)

        # ..CIB (cosmic infrared background)
        Ecib = Ecib * 200e-6 * rm
        NEPcibL = NEP2(det_eff, etendue, Ecib, td3L, Tcib, rm*sc.c)
        NEPcibR = NEP2(det_eff, etendue, Ecib, td3R, Tcib, rm*sc.c)
        NEPcib2 = NEPcibL + NEPcibR
        NEPcib = np.sqrt(NEPcib2)

        # ..from zodiacal light
        Ezodi = Ezodi * np.sqrt(30e-6 * rm)
        NEPzodiL = NEP2(det_eff, etendue, Ezodi, td3L, Tzodi, rm*sc.c)
        NEPzodiR = NEP2(det_eff, etendue, Ezodi, td3R, Tzodi, rm*sc.c)
        NEPzodi2 = NEPzodiL + NEPzodiR
        NEPzodi = np.sqrt(NEPzodi2)

        # ..get total from cosmic radiation
        NEPtotBg2 = NEPcmb2 + NEPcib2 + NEPzodi2
        NEPtotBg = np.sqrt(NEPtotBg2)

        # ..and total from cosmic and instrument backgrounds
        NEPtot2 = NEPtotInst2 + NEPtotBg2
        NEPtot = np.sqrt(NEPtot2)

        self.result['td0L'] = td0L
        self.result['td0R'] = td0R

        self.result['NEPTelR'] = np.sqrt(NEPR1)
        self.result['NEPSecondaryR'] = np.sqrt(NEPR2)
        self.result['NEPBeamsplitterR'] = np.sqrt(NEPR3)
        self.result['NEPStrayR'] = np.sqrt(NEPR4)
        self.result['NEPColdboxR'] = np.sqrt(NEPR5)

        self.result['NEPTelL'] = np.sqrt(NEPL1)
        self.result['NEPSecondaryL'] = np.sqrt(NEPL2)
        self.result['NEPBeamsplitterL'] = np.sqrt(NEPL3)
        self.result['NEPStrayL'] = np.sqrt(NEPL4)
        self.result['NEPColdboxL'] = np.sqrt(NEPL5)

        self.result['NEPtotInst'] = NEPtotInst
        self.result['NEPcmb'] = NEPcmb
        self.result['NEPcib'] = NEPcib
        self.result['NEPzodi'] = NEPzodi
        self.result['NEPtotBg'] = NEPtotBg
        self.result['NEPtot'] = NEPtot

        return self.result

    def __repr__(self):
        return '''
BackgroundNoise:
  NEP Instrument        : {ni}
  NEP Background        : {nb}
  NEP total             : {nt}

  NEP Tel R             : {telR}
  NEP Secondary R       : {secR}
  NEP Beamsplitter R    : {bsR}
  NEP Stray R           : {strayR}
  NEP Cold Box R        : {coldboxR}

  NEP Tel L             : {telL}
  NEP Secondary L       : {secL}
  NEP Beamsplitter L    : {bsL}
  NEP Stray L           : {strayL}
  NEP Cold Box L        : {coldboxL}

  NEP cmb               : {cmb}
  NEP cib               : {cib}
  NEP zodiacal          : {zodi}

'''.format(
          ni=self.result['NEPtotInst'],
          nb=self.result['NEPtotBg'],
          nt=self.result['NEPtot'],

          telR=self.result['NEPTelR'],
          secR=self.result['NEPSecondaryR'],
          bsR=self.result['NEPBeamsplitterR'],
          strayR=self.result['NEPStrayR'],
          coldboxR=self.result['NEPColdboxR'],

          telL=self.result['NEPTelL'],
          secL=self.result['NEPSecondaryL'],
          bsL=self.result['NEPBeamsplitterL'],
          strayL=self.result['NEPStrayL'],
          coldboxL=self.result['NEPColdboxL'],

          cmb=self.result['NEPcmb'],
          cib=self.result['NEPcib'],
          zodi=self.result['NEPzodi'])
