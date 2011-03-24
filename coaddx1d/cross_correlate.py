import numpy as np

#This code was translated from the GHRS IDL library function, CROSS_CORRELATE
#http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library43.html?CROSS_CORRELATE

#HISTORY:
#Version 1  D. Lindler  Sept. 1991

def _cross_correlate(s1, s2, ishift = None, width = None, i1 = None, i2 = None):

    #s1 - first spectrum
    #s2 - second spectrum

    #ishift - approximate offset (default = 0)
    #width - search width (default = 15)
    #i1,i2 - region in first spectrum containing the feature(s)
    #        (default  i1=0, i2=n_elements(s2)-1)

    #outputs:
    #offset: offset of s2 from s1 in data points
    #corr: output correlation vector

    if ishift == None:
        ishift = 0.0
    approx = int((ishift+100000.5) - 100000.0)

    # print("approx: %s" % approx)

    if width == None:
        width = 15

    ns = len(s1)


    if i1 == None:
        i1 = 0
    if i2 == None:
        i2 = ns - 1

    ns2 = ns/2
    width2 = width/2
    it2_start = (i1 - approx + width2) if (i1 - approx + width2) > 0 else 0
    it2_end = (i2 - approx - width2) if (i2 - approx - width2) < (ns - 1) else (ns - 1)
    nt = it2_end - it2_start + 1

    if nt < 1.0:
        raise Exception("cross correlate - region too small, width too large, or ishift too large")

    template2 = s2[it2_start:(it2_end + 1)]

    corr = np.zeros((width))
    mean2 = template2.sum()/nt
    sig2 = np.sqrt(np.sum((template2 - mean2)**2))
    diff2 = template2 - mean2

    for i in range(width):
        it1_start = it2_start - width2 + approx + i
        it1_end = it1_start + nt - 1
        template1 = s1[it1_start:(it1_end + 1)]
        mean1 = template1.sum()/nt
        sig1 = np.sqrt(np.sum((template1 - mean1)**2))
        diff1 = template1 - mean1

        if (sig1 == 0) or (sig2 == 0):
            raise Exception("cross correlate - zero variance computed")

        corr[i] = np.sum(diff1*diff2)/(sig1*sig2)

    maxc = corr.max()
    k = np.where(corr == corr.max())
    k = k[0][0]

    # print(maxc)
    # print(k)

    if (k == 0) or (k == (width - 1.0)):
        raise Exception("cross correlate - maximum on edge of search area")

    kmin = (corr[k - 1] - corr[k])/(corr[k-1] + corr[k+1] - 2.0*corr[k]) - 0.5
    offset = k + kmin - width2 + approx

    return (offset, corr)
