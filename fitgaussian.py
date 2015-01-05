# Gaussian fitting routine adapted from a version on web, written by Andrew Giessel.

import numpy as np
import scipy.optimize

class FitGaussian(object):

    def gaussian(self, height, center_x, center_y, width_x, width_y, rotation):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)
 
        rotation = np.deg2rad(rotation)
        # bug in web version fixed here
        center_xp = center_x * np.cos(rotation) - center_y * np.sin(rotation)
        center_yp = center_x * np.sin(rotation) + center_y * np.cos(rotation)

        def rotgauss(x,y):
            xp = x * np.cos(rotation) - y * np.sin(rotation)
            yp = x * np.sin(rotation) + y * np.cos(rotation)
            g = height*np.exp(-(((center_xp-xp)/width_x)**2+
              ((center_yp-yp)/width_y)**2)/2.)
            return g

        return rotgauss
 
    def moments(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()
        return height, x, y, width_x, width_y, 0.0
 
    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.moments(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) - data)
        p, success = scipy.optimize.leastsq(errorfunction, params)
        return p

