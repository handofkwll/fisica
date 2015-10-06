from __future__ import absolute_import
__author__ = 'rjuanola'
import collections
import numpy as np
import scipy.constants as sc

def black_body(temperature, wavenumber):
    """Function to calculate Planck function.
       temperature - Kelvin
       wavenumber - frequency in cm-1
    """
    freq = wavenumber * sc.c * 100
    bb = 2 * sc.h * np.power(freq, 3) / (np.power(sc.c, 2) / (np.exp((sc.h * freq) / (sc.k * temperature)) - 1.0))
    return bb

def NEP2(eta_spillover, eta_feed, throughput, emissivity, trans2det, temp, freq):
    den = np.exp(sc.h * freq / sc.k / temp) - 1
    a = 1 + (emissivity * trans2det * eta_feed / den)
    b = throughput * emissivity * trans2det * eta_feed * np.power(freq, 4) / den
    cte = 4 * eta_spillover * np.power(sc.h / sc.c, 2)
    d = a * b * (freq[2] - freq[1])
    nep2 = cte * np.sum(d)
    return nep2


class Noise(object):
    """Class to generate the NEP of the simulated observation.
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

        background = self.parameters['substages']['Background']
        coldoptics = self.parameters['substages']['ColdOptics']
        warmoptics = self.parameters['substages']['WarmOptics']

        # Loading background parameters
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

        #Loading cold optics parameters
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

        # Calculating blackbody emission
        B_tel = black_body(Ttel, fts_wn_truncated) * sc.c
        B_stray = black_body(Tstray, fts_wn_truncated) * sc.c
        B_box = black_body(Tbox, fts_wn_truncated) * sc.c
        B_bsm = black_body(Tbsm, fts_wn_truncated) * sc.c
        B_k = black_body(Tk, fts_wn_truncated) * sc.c
        B_sun = black_body(Tsun,fts_wn_truncated) * sc.c

        # Cold optics transmission
        # ..Box left arm to detector 1
        t4L = Weff * NFdicT * np.power(1 - COEmis, NmirrorsL) * DicT * FilT * CoT
        # ..Box right arm to detector 1
        t4R = Weff * NFdicT * np.power(1 - COEmis, NmirrorsR) * DicT * FilT * CoR 

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
        Estray = 0.2 * Etel
        e3 = Estray

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
        det_eff = 0.57
        area = np.pi * np.power(Dtel/2, 2)

        # temporary throughput ? is this correct looks more like the beam but 
        # out by a factor 5 or so
        omegac = np.power(1./rm, 2) / area

        NEPL1 = NEP2(1, det_eff, area * omegac, e1, td1L, Ttel, rm * sc.c)
        NEPL2 = NEP2(1, det_eff, area * omegac, e2, td2L, Ttel, rm * sc.c)
        NEPL3 = NEP2(1, det_eff, area * omegac, ebsm, td12L, Ttel, rm * sc.c)
        NEPL4 = NEP2(1, det_eff, area * omegac, e3, td3L, Tstray, rm * sc.c)
        NEPL5 = NEP2(1, det_eff, area * omegac, Szodi, td0L, Tsun, rm * sc.c)
        NEPL6 = NEP2(1, det_eff, area * omegac, e4L, td4L, Tbox, rm * sc.c)

        NEPR1 = NEP2(1, det_eff, area * omegac, e1, td1R, Ttel, rm * sc.c)
        NEPR2 = NEP2(1, det_eff, area * omegac, e2, td2R, Ttel, rm * sc.c)
        NEPR3 = NEP2(1, det_eff, area * omegac, ebsm, td12R, Tbsm, rm * sc.c)
        NEPR4 = NEP2(1, det_eff, area * omegac, e3, td3R, Tstray, rm * sc.c)
        NEPR5 = NEP2(1, det_eff, area * omegac, Szodi, td0R, Tsun, rm * sc.c)
        NEPR6 = NEP2(1, det_eff, area * omegac, e4R, td4R, Tbox, rm * sc.c)

        NEPL = NEPL1 + NEPL2 + NEPL3 + NEPL4 + NEPL5 + NEPL6
        NEPR = NEPR1 + NEPR2 + NEPR3 + NEPR4 + NEPR5 + NEPR6

        # TOTAL BACKGROUND NEP
        NEPtotInst2 = (np.sum(NEPL) + np.sum(NEPR))
        NEPtotInst = np.sqrt(NEPtotInst2)

        # CMB
        EcmbL = 1 - e1 - e2 - ebsm - e4L
        EcmbR = 1 - e1 - e2 - ebsm - e4R

        NEPcmbL = NEP2(1, det_eff, area*omegac, EcmbL, td3L, Tcmb, rm*sc.c)
        NEPcmbR = NEP2(1, det_eff, area*omegac, EcmbR, td3R, Tcmb, rm*sc.c)
        NEPcmb2 = NEPcmbL + NEPcmbR
        NEPcmb = np.sqrt(NEPcmb2)

        # CIB
        Ecib = Ecib * 200e-6 * rm
        NEPcibL = NEP2(1, det_eff, area*omegac, Ecib, td3L, Tcib, rm*sc.c)
        NEPcibR = NEP2(1, det_eff, area*omegac, Ecib, td3R, Tcib, rm*sc.c)
        NEPcib2 = NEPcibL + NEPcibR
        NEPcib = np.sqrt(NEPcib2)

        #zodi
        Ezodi = Ezodi * np.sqrt(30e-6 * rm)
        NEPzodiL = NEP2(1, det_eff, area*omegac, Ezodi, td3L, Tzodi, rm*sc.c)
        NEPzodiR = NEP2(1, det_eff, area*omegac, Ezodi, td3R, Tzodi, rm*sc.c)
        NEPzodi2 = NEPzodiL + NEPzodiR
        NEPzodi = np.sqrt(NEPzodi2)

        NEPtotBg2 = NEPcmb2 + NEPcib2 + NEPzodi2
        NEPtotBg = np.sqrt(NEPtotBg2)

        NEPtot2 = NEPtotInst2 + NEPtotBg2
        NEPtot = np.sqrt(NEPtot2)

        self.result['NEPtotInst'] = NEPtotInst
        self.result['NEPcmb'] = NEPcmb
        self.result['NEPcib'] = NEPcib
        self.result['NEPzodi'] = NEPzodi
        self.result['NEPtotBg'] = NEPtotBg
        self.result['NEPtot'] = NEPtot

        return self.result


    def __repr__(self):
        return '''
Noise:
  NEP Instrument        : {ni}
  NEP Background        : {nb}
  NEP total             : {nt}
'''.format(
          ni=self.result['NEPtotInst'],
          nb=self.result['NEPtotBg'],
          nt=self.result['NEPtot'])
