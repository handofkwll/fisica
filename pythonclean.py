import fitgaussian
import matplotlib.pyplot as plt
import numpy
import os.path
import pickle
import pythonclean

def hogbom(dirty,
           psf,
           cleanbeam,
           window,
           gain,
           thresh,
           niter):
    """
    Hogbom clean. An algorithm to deconvolve (clean) a messy psf
    from a dirty map.

    Taken from an original written by Bojan Nikolic
    <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk> published
    on the web and licensed under GNU GPL V2 
    (http://www.gnu.org/licenses/gpl.html)

    :param dirty: The dirty image, i.e., the image to be deconvolved

    :param psf: The point spread-function

    :param window: Regions where clean components are allowed. If
    True, think all of the dirty image is assumed to be allowed for
    clean components

    :param gain: The "loop gain", i.e., the fraction of the brightest
    pixel that is removed in each iteration

    :param thresh: Cleaning stops when the maximum of the absolute
    deviation of the residual is less than this value

    :param niter: Maximum number of components to make if the
    threshold "thresh" is not hit
    """

    def overlapIndices(a1, a2, 
                       shiftx, shifty):
        """Utility routine used by hogbom to establish the
        overlap between the psf (a2) centred at shiftx, shifty
        and dirty image (a1). a1 and a2 are assumed to be the
        same shape.
        """

        if shiftx >= 0:
            a1xbeg = shiftx
            a2xbeg = 0
            a1xend = a1.shape[0]
            a2xend = a1.shape[0]-shiftx
        else:
            a1xbeg = 0
            a2xbeg =- shiftx
            a1xend = a1.shape[0]+shiftx
            a2xend = a1.shape[0]

        if shifty >= 0:
            a1ybeg = shifty
            a2ybeg = 0
            a1yend = a1.shape[1]
            a2yend = a1.shape[1]-shifty
        else:
            a1ybeg = 0
            a2ybeg =- shifty
            a1yend = a1.shape[1]+shifty
            a2yend = a1.shape[1]

        return (a1xbeg, a1xend, a1ybeg, a1yend), (a2xbeg, a2xend, a2ybeg, a2yend)

    # arrays to hold clean components, residual and clean image
    comps = numpy.zeros(dirty.shape)
    res = numpy.array(dirty)
    image = numpy.zeros(dirty.shape)

    # if window not set explicitly then clean the whole area
    if window is True:
        window = [slice(0, dirty.shape[0]-1), slice(0, dirty.shape[1]-1)]

    grid = numpy.indices(dirty.shape)

    # loop, looking for clean components
    for i in range(niter):
        # locate brightest pixel in current residual
        argmax = numpy.unravel_index(numpy.fabs(res[window]).argmax(),
          res[window].shape)
        mx = grid[0][window][argmax]
        my = grid[1][window][argmax]

        # add this component to the component list, subtract from the 
        # residual gain*maxval convolved with psf
        mval = res[mx, my]*gain
        comps[mx, my] += mval
        a1o, a2o = overlapIndices(dirty, psf,
                                  mx-dirty.shape[0]/2,
                                  my-dirty.shape[1]/2)
        res[a1o[0]:a1o[1],a1o[2]:a1o[3]] -= psf[a2o[0]:a2o[1],a2o[2]:a2o[3]] * mval
        # also convolve compenent with the clean beam
        image[a1o[0]:a1o[1],a1o[2]:a1o[3]] += cleanbeam[a2o[0]:a2o[1],a2o[2]:a2o[3]] * mval

        if numpy.fabs(res).max() < thresh:
            break

    return comps, res, image

def smoothness(comps):
    """NOT FOR USE. This routine is part of the implementation of a scheme suggested by
    T.J.Cornwell 1983, A&A, 121, 281 'A method for stabilizing the clean algorithm'.
    It is not working yet, and may not be necessary for reasonable uv coverage.
    """
 
    s = -numpy.sum(comps**2)
    return s

def modified(comps, dirty, sigma):
    """NOT FOR USE. This routine is part of the implementation of a scheme suggested by
    T.J.Cornwell 1983, A&A, 121, 281 'A method for stabilizing the clean algorithm'.
    It is not working yet, and may not be necessary for reasonable uv coverage.
    """

    imshape = numpy.shape(comps)
    irange = numpy.arange(imshape[1])
    jrange = numpy.arange(imshape[0])

    h = smoothness(comps)

    # h <= 0
    dik = numpy.sqrt(numpy.mean(-h))
    if dik==0:
        dik = 1.0

    dh_by_dik = numpy.zeros(imshape)

    for i in irange:
        for j in jrange:
            dh = h + comps[j,i]**2 - (comps[j,i] + dik)**2
            dh_by_dik[j,i] = dh / dik 

    alpha_opt = sigma**2 / (2 * numpy.mean(comps * dh_by_dik))

    d_modified = dirty + alpha_opt * dh_by_dik

    return d_modified

def cleantest(psffile='dirtybeam_100uvpoints.pickle'):
    """Procedure to test hogbom. The dirty beam is read in from
    psffile, and convolved with a simple pattern to make a dirty
    map.

    'hogbom' is then called to clean away the psf.

    Various results are then plotted to file 'cleantest.png'.
    """

    # obtain the psf
    filename = os.path.join(os.path.dirname(__file__), 'data',
      psffile)
    print 'Reading psf from: %s' % filename
    f = open(filename)
    psf = pickle.load(f)
    f.close()

    imshape = numpy.shape(psf)

    # fit Gaussian to centre of psf and use this as the 'clean beam'
    fitter = fitgaussian.FitGaussian()
    psf_centre = numpy.array(psf[imshape[0]/2-5:imshape[0]/2+5,
      imshape[1]/2-5:imshape[1]/2+5])
    p = fitter.fitgaussian(psf_centre)

    # construct the clean beam
    cp = (1.0, float(imshape[0])/2.0, float(imshape[1])/2.0, p[3], p[4], p[5])
    rotgauss = fitter.gaussian(*cp)
    cleanbeam = numpy.fromfunction(rotgauss, imshape)

    # construct a dirty image from a simple ideal image convolved with 
    # the psf.
    dirty = numpy.zeros(imshape)
    jcen = imshape[1]/2
    icen = imshape[0]/2

    for j in range(imshape[1]):
        for i in range(imshape[0]):
            
            if ((i-icen)**2 + (j-jcen)**2) < 400:
                pixel = (i, j)  
# debug                pixel = (imshape[0]/2, imshape[1]/2)  
# debug                pixel = (5,5)

                d_imin = pixel[0] - imshape[0]/2
                d_imax = pixel[0] + imshape[0]/2 - 1
                d_jmin = pixel[1] - imshape[1]/2
                d_jmax = pixel[1] + imshape[1]/2 - 1

                d_imin = max(0, d_imin)
                d_imax = min(imshape[0]-1, d_imax)
                d_jmin = max(0, d_jmin)
                d_jmax = min(imshape[1]-1, d_jmax)

                p_imin = imshape[0]/2 - pixel[0]
                p_imax = imshape[0] - 1 + imshape[0]/2 - pixel[0]
                p_jmin = imshape[1]/2 - pixel[1]
                p_jmax = imshape[1] - 1 + imshape[1]/2 - pixel[1]

                p_imin = max(0, p_imin)
                p_imax = min(imshape[0]-1, p_imax)
                p_jmin = max(0, p_jmin)
                p_jmax = min(imshape[1]-1, p_jmax)

                dirty[d_imin:d_imax, d_jmin:d_jmax] += psf[p_imin:p_imax, p_jmin:p_jmax]

    for cln in range(1):
        print 'test clean loop', cln
        print '..cleaning'
        comps, res, image = pythonclean.hogbom(
                                   dirty=dirty,
                                   psf=psf,
                                   cleanbeam=cleanbeam,
                                   window=[slice(imshape[0]/4, imshape[0]*3/4),
                                   slice(imshape[1]/4, imshape[1]*3/4)],
                                   gain=0.1,
                                   thresh=0.01,
                                   niter=1000)

# following commented out lines were experimenting with one method
# developed to minimise 'ringing' in the cleaned map. This method
# was not got working, and may not be needed if uv coverage is good
# enough.
#        print '..new dirty image'
#        # modify the dirty image
#        dirty = modified(comps, dirty, sigma=0.01)

    # plot the test results.
    plt.figure()

    plt.subplot(221)
    plt.imshow(psf, interpolation='nearest', origin='lower',
      aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
      vmax=numpy.max(psf), vmin=numpy.min(psf))
    plt.colorbar(orientation='vertical')
    plt.axis('image')
    plt.title('psf')

    plt.subplot(222)
    plt.imshow(dirty, interpolation='nearest', origin='lower',
      aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
      vmax=numpy.max(dirty), vmin=numpy.min(dirty))
    plt.colorbar(orientation='vertical')
    plt.axis('image')
    plt.title('dirty')

    plt.subplot(223)
    plt.imshow(image, interpolation='nearest', origin='lower',
      aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
      vmax=numpy.max(image), vmin=numpy.min(image))
    plt.colorbar(orientation='vertical')
    plt.axis('image')
    plt.title('image')

    plt.subplot(224)
    plt.imshow(res, interpolation='nearest', origin='lower',
      aspect='equal', extent=[0, imshape[0], 0, imshape[1]],
      vmax=numpy.max(res[imshape[0]/4:imshape[0]*3/4,
      imshape[1]/4:imshape[1]*3/4]),
      vmin=numpy.min(res[imshape[0]/4:imshape[0]*3/4,
      imshape[1]/4:imshape[1]*3/4]))
    plt.colorbar(orientation='vertical')
    plt.axis('image')
    plt.title('residual')

    filename = 'cleantest.png'
    plt.savefig(filename)
    plt.close()
