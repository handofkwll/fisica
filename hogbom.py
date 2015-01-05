import numpy 

def image(himage, psf, dirty, niter, thresh):
    himage.open(psf)
    psfdata = himage.getchunk()
    print 'psf read, max val %s' % numpy.max(psfdata)
    psfdata /= numpy.max(psfdata)
    print 'psf normalised'

#    if not dirty: 
#        sky = numpy.zeros(numpy.shape(psf))
#        dirty = numpy.zeros(numpy.shape(psf))
    
#        sky[70,70]=5000
#        for i in range(numpy.shape(psf)[0]):
#            for j in range(numpy.shape(psf)[1]):
#                if pow((i-20),2) + pow((j-20),2) < 100:
#                    sky[i,j] = 10000

#        for i in range(numpy.shape(psf)[0]):
#            for j in range(numpy.shape(psf)[1]):
#                a1o, a2o=overlapIndices(sky, psf,
#                                    i+1-sky.shape[0]/2,
#                                    j+1-sky.shape[1]/2)
#                dirty[a1o[0]:a1o[1],a1o[2]:a1o[3]]+=psf[a2o[0]:a2o[1],a2o[2]:a2o[3]]*sky[i,j]
#    else:
#        himage.open(dirty)
#        dirtydata = himage.getchunk()
#    himage.fromarray(outfile='dirty.im', pixels=dirty, overwrite=True) 
#    print 'dirty read'

    himage.open(dirty)
    dirtydata = himage.getchunk()
    print 'dirty read'

    comps, res = hogbom(dirtydata, psfdata, True, 0.05, thresh, niter)
    print 'clean completed'

    himage.fromarray(outfile='residual.im', pixels=res, overwrite=True) 
    himage.fromarray(outfile='comps.im', pixels=comps, overwrite=True) 
    himage.fromarray(outfile='image.im', pixels=comps+res, overwrite=True) 
    himage.close()


# Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk> 
# Initial version August 2010
#
# This file is part of pydeconv. This work is licensed under GNU GPL
# V2 (http://www.gnu.org/licenses/gpl.html)
"""
Clean based deconvolution, using numpy
"""

def overlapIndices(a1, a2, 
                   shiftx, shifty):
    if shiftx >=0:
        a1xbeg=shiftx
        a2xbeg=0
        a1xend=a1.shape[0]
        a2xend=a1.shape[0]-shiftx
    else:
        a1xbeg=0
        a2xbeg=-shiftx
        a1xend=a1.shape[0]+shiftx
        a2xend=a1.shape[0]

    if shifty >=0:
        a1ybeg=shifty
        a2ybeg=0
        a1yend=a1.shape[1]
        a2yend=a1.shape[1]-shifty
    else:
        a1ybeg=0
        a2ybeg=-shifty
        a1yend=a1.shape[1]+shifty
        a2yend=a1.shape[1]

    return (a1xbeg, a1xend, a1ybeg, a1yend), (a2xbeg, a2xend, a2ybeg, a2yend)

        

def hogbom(dirty,
           psf,
           window,
           gain,
           thresh,
           niter):
    """
    Hogbom clean

    :param dirty: The dirty image, i.e., the image to be deconvolved

    :param psf: The point spread-function

    :param window: Regions where clean components are allowed. If
    True, thank all of the dirty image is assumed to be allowed for
    clean components

    :param gain: The "loop gain", i.e., the fraction of the brightest
    pixel that is removed in each iteration

    :param thresh: Cleaning stops when the maximum of the absolute
    deviation of the residual is less than this value

    :param niter: Maximum number of components to make if the
    threshold "thresh" is not hit
    """
    comps=numpy.zeros(dirty.shape)
    res=numpy.array(dirty)
    if window is True:
        window=numpy.ones(dirty.shape,
                          numpy.bool)
    for i in range(niter):
        mx, my=numpy.unravel_index(numpy.fabs(res[window]).argmax(), res.shape)
        print 'mx, my', mx, my
        mval=res[mx, my]*gain
        print 'mval', mval
        comps[mx, my]+=mval
        a1o, a2o=overlapIndices(dirty, psf,
                                mx+1-dirty.shape[0]/2,
                                my+1-dirty.shape[1]/2)
#                                mx-dirty.shape[0]/2,
#                                my-dirty.shape[1]/2)
        print 'a1o', a1o
        print 'a2o', a2o
        res[a1o[0]:a1o[1],a1o[2]:a1o[3]]-=psf[a2o[0]:a2o[1],a2o[2]:a2o[3]]*mval
        print '59,56', res[59,56]
        if numpy.fabs(res).max() < thresh:
            break
    return comps, res
