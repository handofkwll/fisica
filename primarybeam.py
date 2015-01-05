from __future__ import absolute_import

import collections
import numpy as np
import pp

import matplotlib.pyplot as plt


class PrimaryBeam(object):
    """Primary beam simulator.
    """
    def __init__(self, job_server=None, beamfile='Sky_Beam_Example.dat'):
        self.beamfile = beamfile
        self.result = collections.OrderedDict()

        if job_server is None:
            # start parallel python (pp), find the number of CPUS available
            self.job_server = pp.Server(ppservers=())
            print 'PrimaryBeam starting pp with %s workers' % \
              self.job_server.get_ncpus()
        else:
            self.job_server = job_server

    def readfile(self):
        f = open(self.beamfile, 'r')

        imsize = f.readline()
        imsize = imsize.split()
        imsize = [int(dim) for dim in imsize]

        # could check headings eventually
        headings = f.readline().split()

        u = set()
        v = set()
        xamp_store = {}
        xphase_store = {}
        yamp_store = {}
        yphase_store = {}
        zamp_store = {}
        zphase_store = {}

        words = f.readline().split()
        while words:
            numbers = map(float, words)
            assert len(numbers)==8, 'unexpected number of numbers on line: %s' % line

            u.update([numbers[6]])
            v.update([numbers[7]])

            xamp_store[(numbers[6], numbers[7])] = numbers[0]
            xphase_store[(numbers[6], numbers[7])] = numbers[1]
            yamp_store[(numbers[6], numbers[7])] = numbers[2]
            yphase_store[(numbers[6], numbers[7])] = numbers[3]
            zamp_store[(numbers[6], numbers[7])] = numbers[4]
            zphase_store[(numbers[6], numbers[7])] = numbers[5]

            words = f.readline().split()

        f.close()

        u = list(u)
        u.sort()
        self.u = np.array(u)
        v = list(v)
        v.sort()
        self.v = np.array(v)

        xamp = np.zeros([len(u), len(v)])
        xphase = np.zeros([len(u), len(v)])
        yamp = np.zeros([len(u), len(v)])
        yphase = np.zeros([len(u), len(v)])
        zamp = np.zeros([len(u), len(v)])
        zphase = np.zeros([len(u), len(v)])
        for k in xamp_store.keys():
            xamp[self.v==k[1], self.u==k[0]] = xamp_store[k]
            xphase[self.v==k[1], self.u==k[0]] = xphase_store[k]
            yamp[self.v==k[1], self.u==k[0]] = yamp_store[k]
            yphase[self.v==k[1], self.u==k[0]] = yphase_store[k]
            zamp[self.v==k[1], self.u==k[0]] = zamp_store[k]
            zphase[self.v==k[1], self.u==k[0]] = zphase_store[k]

        self.xfield = np.zeros([len(u), len(v)], np.complex)
        self.xfield.real = xamp * np.cos(xphase)
        self.xfield.imag = xamp * np.sin(xphase)
        self.yfield = np.zeros([len(u), len(v)], np.complex)
        self.yfield.real = yamp * np.cos(yphase)
        self.yfield.imag = yamp * np.sin(yphase)
        self.zfield = np.zeros([len(u), len(v)], np.complex)
        self.zfield.real = zamp * np.cos(zphase)
        self.zfield.imag = zamp * np.sin(zphase)

    def plots(self):
        plt.figure()

        plt.subplot(231)
        plt.imshow(np.abs(self.xfield), interpolation='nearest', origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.abs(self.xfield))*1.1,
          vmin=np.min(np.abs(self.xfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test X amp')

        plt.subplot(232)
        plt.imshow(np.angle(self.xfield), interpolation='nearest', origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.angle(self.xfield))*1.1,
          vmin=np.min(np.angle(self.xfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test X phase')

        plt.subplot(233)
        plt.imshow(np.abs(self.yfield), interpolation='nearest', origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.abs(self.yfield))*1.1,
          vmin=np.min(np.abs(self.yfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test Y amp')

        plt.subplot(234)
        plt.imshow(np.angle(self.yfield), interpolation='nearest', origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.angle(self.yfield))*1.1,
          vmin=np.min(np.angle(self.yfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test Y phase')

        plt.subplot(235)
        plt.imshow(np.abs(self.zfield), interpolation='nearest', origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.abs(self.zfield))*1.1,
          vmin=np.min(np.abs(self.zfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test Z amp')

        plt.subplot(236)
        plt.imshow(np.angle(self.zfield), interpolation='nearest',
          origin='lower',
          aspect='equal', extent=[self.u[0], self.u[-1],
          self.v[0], self.v[-1]], vmax=np.max(np.angle(self.zfield))*1.1,
          vmin=np.min(np.angle(self.zfield))*0.9)
        plt.colorbar(orientation='vertical')
        plt.axis('image')
        plt.title('Test Z phase')

        filename = 'primary.png'
#        filename = sanitize(filename)
#        filename = os.path.join(context['dirname'], filename)
        plt.savefig(filename)
        plt.close()

    def __repr__(self):
        return 'FISICA PrimaryBeam object'

