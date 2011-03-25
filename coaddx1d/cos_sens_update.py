from numpy import array, ndarray

# Documentation from original program:

# pass in wave, segment, and grating and get out the smov sensitivity
# update USE the lv option for data taken after Aug 13, 2009 divide
# the returned sensvector value times your SMOV flux to get updated
# values

# Original IDL code written by Steve Penton.

def _cos_sens_update(wave, flux, segment = "FUVA", lv = False, grating = "G130M", p = None, sensvector = None):

    if (type(wave) and type(flux)) != ndarray:
        raise Exception("wave and flux inputs must be of type numpy.array")

    if (type(segment) and type(grating)) != str:
        raise Exception("segment and grating inputs must be strings")

    if type(lv) != bool:
        raise Exception("lv must be of boolean type")

    nw = len(wave)
    nf = len(flux)

    if nw != nf:
        raise Exception("Error: the wavelength and flux vectors must have the same size")

    useg = segment.upper().strip()
    ugrat = grating.upper().strip()

    if ugrat == "G130M":
        pa=[2736.3667,-6.5024309,0.0040286034,1.2969438e-06,-2.0995578e-09,5.5387341e-13]
        pb=[7385.5333,-13.271672,-0.0048998948,2.4903635e-05,-1.8383126e-08,4.2931604e-12]
        pa_lv=[913.33789,5.3526649,-0.020921567,2.5117444e-05,-1.2848511e-08,2.4239719e-12]
        pb_lv=[6566.5360,-12.140919,-0.0032876628,2.0892240e-05,-1.5698942e-08,3.6928369e-12]
    elif ugrat == "G160M":
        pa=[-3279.1462,5.9503510,-0.0025139167,-1.2436785e-06,1.1917261e-09,-2.4044604e-13]
        pb=[-5080.8995,9.1067530,-0.0015959567,-6.0123991e-06,4.3706627e-09,-8.9840028e-13]
        pa_lv=[-1378.5339,2.6580205,-0.0014300817,-1.8917433e-07,3.6880377e-10,-8.1514001e-14]
        pb_lv=[-7073.9153,12.539628,-0.0019500167,-8.5522165e-06,6.1290053e-09,-1.2525614e-12]
    elif ugrat == "G140L":
        raise Exception("Error: G140L data are not yet supported.")
    else:
        raise Exception("Error: grating %s is incorrect" % ugrat)

    if useg == "FUVA":
        if lv == True:
            p = pa_lv
        else:
            p = pa
    else:
        if lv == True:
            p = pb_lv
        else:
            p = pb

    sensvector = array([0.0]*nw)
    
    for n, i in enumerate(p):
        sensvector += float(i)*(wave**n)

    return flux/sensvector
