"""
This module provides the coadd function.

HISTORY: written by Charles Danforth (danforth@casa.colorado.edu)
substantial input and testing by Brian Keeney, Kevin France, Yangsen
Yao, Hao Yang, Steve Penton (U. Colorado) and Anand Narayanan
(Wisconsin)

The IDL source code is located here:
http://casa.colorado.edu/~danforth/science/cos/coadd_x1d.pro

Translated into python by Nico Nell (nicholas.nell@colorado.edu)
"""


import os, fnmatch
import numpy as np
import pyfits
import idlsave
import csv
import pickle
from datetime import datetime
from scipy.signal import medfilt
from cos_sens_update import _cos_sens_update
from cross_correlate import _cross_correlate


# TODO We need the wireflat_iter.idl file for the newflat option...
# TODO: test this and previous methods to check exposure time. This
# method (vs the concatenates) has some peaks that are probably caused
# by /0

# Ignore divide by 0 math errors.
np.seterr(divide='ignore', invalid='ignore')

class COSSpectrum(object):
    '''
    HST-COS Spectrum object.
    
    It serves mainly as a convenient container for the output of coaddition.
    '''

    def __init__(self, readme, targname, fitshdr, files, grating, side,
                 cenwave, fppos, dateobs, exptin, wavein, fluxin, errin,
                 wave, flux, err, exptime):

        self.readme = readme
        self.targname = targname
        self.fitshdr = fitshdr
        self.files = files
        self.grating = grating
        self.side = side
        self.cenwave = cenwave
        self.fppos = fppos
        self.dateobs = dateobs
        self.exptin = exptin
        self.wavein = wavein
        self.fluxin = fluxin
        self.errin = errin
        self.wave = wave
        self.flux = flux
        self.err = err
        self.exptime = exptime
        

def coadd(files=None, path='.', chan=1, method=1,
          no_flat=False, no_flange=False, update_flux=False,
          no_align=False, scale=None, indivshift=None,
          ubershift=None, savefile=None, saveformat="dat",
          newflat=None, verbose=True):

    '''
    Coadd Cosmic Origins Spectrograph x1d files and return a COSSpectrum object

    Keyword Arguments: 
    
        *files*: [ list of strings ]
        List of files to coadd. If this option is not specified then coadd will
        use all x1d files in the given path.

        *path*: [ str ]
        Path where the x1d files to coadd are.
        Default: current working directory.

        *chan*: [ int ]
        This number corresponds to the gratings used:
        
        1. G130M (Default)
        2. G160M
        3. Both M gratings (G130M + G160M)
        4. G140L

        *method*: [ int ]
        This number corresponds to the method used to coadd:
        
        - -1 - simple mean
        - 0 -
        - 1 - modified exposure weighting (Default)
        - 2 - error squared weighting
        - 3 - signal-to-noise squared weighting

        *no_flat*: [ bool ]
        Do not perform first-order flat fielding (aka gridwire removal)
        Default is false

        *no_flange*: [ bool ]
        Do not de-weight pixels at the edges of the detectors in each exposure
        Default is false

        *update_flux*: [ bool ]
        Correct early flux calibration. This should no longer be necessary
        Default is false

        *no_align*: [ bool ]
        Do not cross-correlate exposures
        Default is false

        *scale*: [ list of floats ]
        Set of scaling factors to multiply input fluxen, errors by. Should be one
        scale per input file.
        Default is None

        *indivshift*: [ list of floats ]
        Exposure/side specific shift (in angstroms). Vector must be [n,2] where n
        is the number of files.

        *ubershift*: [ float ]
        Constant shift (in angstroms) to shift the wavelength by.

        *savefile*: [ str or bool ]
        Name for a file to save the coadded data to. If set to True then a
        default name of "<target>-<grating>" is used.
        (Default is None)

        *saveformat*: [ str ]
        Format to use for the save file:

        - "pickle" - Save the data as a pickled python object.
        - "dat" - Save the file as an ascii space separated file (Default).
        - "csv" - Save the file as a comma separated value file.

        *newflat*: [ bool ]
        NOT YET SUPPORTED!
        Use the itterative flatfield files
        Default is false

        *verbose*: [ bool ]
        This is used to show or hide stdout messages from this program.
        Default is True

    Returns:
        
        *spectrum*: [ COSSpectrum ]
        a COSSpectrum object containing all the necessary information.
    
    '''
    
    # wirefile list: (kept with the python module)
    wirefile = ["wireflat_a.idl", "wireflat_b.idl"]

    # name used on plots
    channame=['G130M', 'G160M', 'G130M + G160M', 'G140L']
    # detector names
    segname=['FUVA','FUVB']
    # grating dispersion (Angstrom/pix)
    disp=[0.00997, 0.01223, -1.0, 0.0803]

    minxcorwave = [[1330,1664,-1,1300],[1255,1520,-1,1000]]
    maxxcorwave = [[1340,1676,-1,1370],[1266,1533,-1,1100]]

    # error-flanging scale height, pixels
    flange_size = 100.0
    # local S/N threshhold for dynamic weighting
    # This is never used... ?
    # snthresh = 5.0

    if method > 3 or method < -1:
        raise Exception("keyword: method must be between -1 and 3")

    if chan == 3:
        #G130M + G160M
        chanind = [0,1]
    else:
        chanind = int(chan - 1)

    # Add a / to the end of path if the user didn't add it..
    if path[-1] != os.path.sep:
        path = path + os.path.sep

    if files != None:
        nfiles = len(files)
    else:
        #this will become an array.
        files = []
        for filename in os.listdir(path):
            if fnmatch.fnmatch(filename, '*_x1d.fits'):
                files.append(filename)

        if len(files) == 0:
            raise Exception("No x1d files found in path %s" % os.path.abspath(path))
        else:
            nfiles = len(files)

        files = np.array(files)

    if scale == None:
        scale = [1.0]*len(files)

    grating = np.array(['']*nfiles, dtype='S50')
    cenwave = np.array([0.0]*nfiles)
    fppos = np.array([0.0]*nfiles)
    filenote = np.array(['']*nfiles, dtype='S50')
    
    for n, i in enumerate(files):
        hdulist = pyfits.open("".join((path,i)))
        grating[n] = hdulist[0].header['OPT_ELEM']
        cenwave[n] = hdulist[0].header['CENWAVE']
        fppos[n] = hdulist[0].header['FPPOS']
        hdulist.close()

    if type(chanind) == int:
        g = np.where(grating == channame[chanind])
        g = g[0]
    else:
        g = set(np.where(grating == channame[chanind[0]])[0]).union(
            set(np.where(grating == channame[chanind[1]])[0]))
        g = np.array(list(g))

    if len(g) > 0:
        grating = grating[g]
        cenwave = cenwave[g]
        fppos = fppos[g]
        files = files[g]
        nfiles = len(g)
    else:
        raise Exception("No files match requested grating: %s" % channame[chan - 1])

    spec = pyfits.open("".join((path,files[0])))
    fitshdr = spec[1].header
    targname = spec[0].header['TARGNAME']
    # index #3 is wavelength
    nwave = len(spec[1].data[0][3])
    spec.close()
    
    boxsize = 100
    zeroes = np.zeros((nwave))
    wavein = np.zeros((nwave, nfiles, 2))
    fluxin = np.zeros((nwave, nfiles, 2))
    errin = np.zeros((nwave, nfiles, 2))
    pix2waveshift = np.zeros((nfiles, 2))
    exptin = np.zeros((nwave, nfiles, 2))
    dateobs = np.array([None]*nfiles, dtype='S50')
    timeobs = np.array([None]*nfiles, dtype='S50')
    
    for n, i in enumerate(files):
        spec = pyfits.open("".join((path,i)))
        dateobs[n] = spec[1].header['DATE-OBS']
        timeobs[n] = spec[1].header['TIME-OBS']
        pix2waveshift[n,:] = -1.0*np.array([spec[1].header['SHIFT1A'], spec[1].header['SHIFT1B']])

        #parse segments A, B into wave, flux, err vectors
        for j in range(2):
            for k in range(nwave):
                lo = 0 if 0 > (k - boxsize/2) else (k - boxsize/2)
                hi = (nwave - 1) if (nwave - 1) < (k + boxsize/2) else (k + boxsize/2)
                aa = np.where(spec[1].data[j][4][lo:hi + 1] == 0)
                naa = len(aa[0])
                zeroes[k] = naa

            aa = np.where(zeroes == 0.0)
            naa = len(aa[0])

            waveintmp = spec[1].data[j][3][aa[0][0]:aa[0][naa - 1] + 1]

            # UPDATE FLUX CALIBRATION ,cos_sens_update can't (yet)
            # handle G140L data, Data calibrated under older versions
            # of CalCOS contained dramatically-incorrect flux
            # calibration as opposed to the only slightly-incorrect
            # flux calibration now. This step shouldn't be necessary
            # for more current reductions.

            if update_flux == True:
                fluxintmp = _cos_sens_update(spec[1].data[j][3], spec[1].data[j][4],
                                             segment = segname[j], grating = grating[n].tostring())
                errintmp =  _cos_sens_update(spec[1].data[j][3], spec[1].data[j][5],
                                             segment = segname[j], grating = grating[n].tostring())
            else:
                fluxintmp = spec[1].data[j][4]
                errintmp = spec[1].data[j][5]

            exptintmp = np.array([spec[1].header['EXPTIME']]*len(fluxintmp))

            # limit errors - in error-weighted coaddition, undue
            # weight is given to points with ridiculously small
            # errors.  For data where there are significant regions of
            # F~0, use the methods=+-1 keyword.

            minerr = 1.0e-17
            errin[np.where(errin < minerr)] = minerr
            fluxintmp = fluxintmp*scale[n]
            errintmp = errintmp*scale[n]

            if no_flat == False:
                # NEW_FLAT
                # NEED TO GET WIREFLAT_ITER.idl from Danforth.

                # Add new flat support...
                # grab the wirefile...
                wifile = idlsave.read("".join((os.path.split(
                    os.path.abspath(__file__))[0] + os.path.sep, wirefile[j])), verbose=False)
                    #os.path.abspath(coaddx1d.__path__) +
                    #                           "/", wirefile[j])))
                flat = np.roll(wifile['pseudoflat'], int(np.round(pix2waveshift[n,j])))
                fluxintmp = fluxintmp/flat
                fatflat = flat*0
                errbox = 7 #pixels, half-width

                for k in range(int(errbox), int(len(flat) - errbox)):
                    fatflat[k] = flat[k - errbox:k + errbox + 1].min()

                errintmp = errintmp/(fatflat**2)
                exptintmp = exptintmp*(fatflat**2)

            fluxintmp = fluxintmp[aa[0][0]:aa[0][naa - 1] + 1]
            errintmp = errintmp[aa[0][0]:aa[0][naa - 1] + 1]
            exptintmp = exptintmp[aa[0][0]:aa[0][naa - 1] +1]

            if no_flange == False:
                npix = len(waveintmp)
                # This creates an array of INTS ... should be float
                # edgedist = np.arange(npix)
                # fixed.
                edgedist = np.array(range(npix), dtype=float)

                flange = 1.0 + np.exp(-1.0*(edgedist + 1.0)/flange_size)
                flange[np.where(flange < (np.exp(-1.0*(npix - 2.0 - edgedist)/flange_size)))] = \
                                       (np.exp(-1.0*(npix - 2.0 - edgedist)/flange_size))

                errintmp = errintmp*flange
                exptintmp = exptintmp/(flange**2)

            # copy temporary arrays to official business.
            wavein[0:len(waveintmp), n, j] = waveintmp
            fluxin[0:len(waveintmp), n, j] = fluxintmp
            errin[0:len(waveintmp), n, j] = errintmp
            exptin[0:len(waveintmp), n, j] = exptintmp

            # pad out the rest of the wavelength range to the maximum
            # wavelength+1A (prevents wave=0)
            wavein[len(waveintmp):, n, j] = waveintmp.max() + 0.01
            fluxin[len(waveintmp):, n, j] = 0.0
            errin[len(waveintmp):, n, j] = 0.0

            # filter out single-pixel down-spikes in error by
            # comparing values to seven-pixel medians.  This takes
            # care of some of it, but there are still problems. 2/2/10

            smootherrin = medfilt(errin[:, n , j], 7)

            for k in range(3,len(waveintmp) - 4):
                if errin[k, n, j] <= 0.1*smootherrin[k]:
                    errin[k, n, j] = np.median(errin[k-3:k+4, n, j])
            
        spec.close()

    xshift = np.zeros((nfiles, 2))

    if indivshift == True:
        no_align = True

    if nfiles > 1 and no_align == False:
        bestgrate = np.array([1309.0, 1600.0, -1.0, 1230.0])
        if type(chanind) == int:
            chanind_dum = [chanind]
        else:
            chanind_dum = chanind

        for i in chanind_dum:
            # ref = np.concatenate((np.where(grating == channame[i]),
            #                       np.where(cenwave == bestgrate[i]),
            #                       np.where(fppos == 3)),axis = 1)
            ref = np.intersect1d(np.where(grating == channame[i])[0],
                                 np.where(cenwave == bestgrate[i])[0])
            ref = np.intersect1d(ref, np.where(fppos == 3)[0])
            
            #ref = (grating == channame[i]) & (cenwave == bestgrate[i]) & (fppos == 3)
            #ref = ref[0]

            if len(ref) == 0:
                ref = np.intersect1d(np.where(grating == channame[i])[0],
                                     np.where(fppos == 3)[0])
                #ref = (grating == channame[i]) & (fppos == 3)
                # ref = np.concatenate((np.where(grating == channame[i]),
                #                       np.where(fppos == 3)), axis = 1)
                #ref = ref[0]
            if len(ref) == 0:
                ref = np.where(grating == channame[i])
                ref = ref[0]
            if len(ref) > 1:
                if verbose:
                    print(" --- ")
                    print("NOTE: multiple reference exposures available for grating: %s" %
                          grating[ref[0]])

                if verbose:
                    for j in ref:
                        print(str(j).strip() + "-" + " ".join((files[j],str(grating[j]),
                                                               str(cenwave[j]),str(fppos[j]),
                                                               str(exptin[:,j,0].max()),
                                                               str(np.median(fluxin[:,j,0])),
                                                               str(exptin[:,j,0].max()),
                                                               str(np.median(fluxin[:,j,1])),
                                                               str(dateobs[j]),
                                                               str(timeobs[j]))))
                        
                if verbose:
                    print("Using exposure 0 above as alignment reference")
                    print("Re-run code with explicit file list if this is non-optimal.")
                    print(" --- ")
                
                
            ref = ref[0]
            filenote[ref] = filenote[ref] + "*"
            thisgrat = np.where(grating == grating[ref])

            xcor_width = 30
            for j in range(2):
                for k in thisgrat[0]:
                    # refrange does not need to be reassigned every
                    # loop. no k dependence...
                    # refrange = (wavein[:, ref, j] >= minxcorwave[i][j]) & \
                    #            (wavein[:, ref, j] <= maxxcorwave[i][j])
                    # comprange = (wavein[:, k, j] >= minxcorwave[i][j]) & \
                    #             (wavein[:, k, j] <= maxxcorwave[i][j])

                    # Note reverse from idl's [column, row] indexing...
                    refrange = (wavein[:, ref, j] >= minxcorwave[j][i]) & \
                               (wavein[:, ref, j] <= maxxcorwave[j][i])
                    comprange = (wavein[:, k, j] >= minxcorwave[j][i]) & \
                                (wavein[:, k, j] <= maxxcorwave[j][i])

                    # The following did not work because IDL uses [column, row] indexing..
                    # refrange = np.concatenate((np.where(wavein[:,ref,j] >= minxcorwave[i][j]),
                    #                            np.where(wavein[:,ref,j] <= maxxcorwave[i][j])),
                    #                           axis = 1)
                    # comprange = np.concatenate((np.where(wavein[:,k,j] >= minxcorwave[i][j]),
                    #                            np.where(wavein[:,k,j] <= maxxcorwave[i][j])),
                    #                           axis = 1)
                    
                    refx = wavein[refrange, ref, j]
                    refy = fluxin[refrange, ref, j]

                    compx = wavein[comprange, k, j]
                    compy = fluxin[comprange, k, j]

                    (refx, refy, compx, compy) = _herczeg(refx,
                                                          refy,
                                                          compx,
                                                          compy,
                                                          verbose = verbose)

                    (shift, corr) = _cross_correlate(refy,
                                                     compy,
                                                     width = xcor_width)

                    xshift[k, j] = shift*disp[i]
                    wavein[:, k, j] = wavein[:, k, j] + xshift[k, j]
    else:
        if indivshift == None:
            xshift = np.zeros((nfiles, 2))
        else:
            xshift = indivshift
            for i in range(nfiles):
                for j in range(1):
                    if verbose:
                        print("Individual wavelength shift for %d, %d, %s, %s" %
                              (i,j,xshift[i,j],indivshift[i,j]))
                    wavein[:, i, j] = wavein[:, i, j] + xshift[i, j]

    #print summary
    #do I want to do this? ... probably...
    if verbose:
        print("COADD_X1D: Exposure Summary")
        for i in range(nfiles):
            print(str(i) + "-" + " ".join((str(files[i]),
                                           str(grating[i]),
                                           str(cenwave[i]),
                                           str(fppos[i]),
                                           format(exptin[:,i,0].max(), ".1e"),
                                           format(np.median(fluxin[:,i,0]), ".2e"),
                                           format(xshift[i,0], ".3f"),
                                           format(exptin[:,i,1].max(), ".1e"),
                                           format(np.median(fluxin[:,i,1]), ".2e"),
                                           format(xshift[i,1], ".3f"),
                                           str(dateobs[i]),
                                           str(timeobs[i]),
                                           str(filenote[i]))))

        print("(* alignment reference exposure)\n")

    # collapse sides A and B into the same vector (remove third
    # dimension from the arrays)

    wavein = np.reshape(wavein, (wavein.shape[0], wavein.shape[1]*wavein.shape[2]))
    fluxin = np.reshape(fluxin, (fluxin.shape[0], fluxin.shape[1]*fluxin.shape[2]))
    errin = np.reshape(errin, (errin.shape[0], errin.shape[1]*errin.shape[2]))
    exptin = np.reshape(exptin, (exptin.shape[0], exptin.shape[1]*exptin.shape[2]))
    grating = np.array([grating, grating]).reshape(grating.shape[0]*2)
    cenwave = np.array([cenwave, cenwave]).reshape(cenwave.shape[0]*2)
    fppos = np.array([fppos, fppos]).reshape(fppos.shape[0]*2)
    side = np.array([['A']*nfiles, ['B']*nfiles])
    side.reshape(nfiles*2)
    nfiles = nfiles*2

    # Define array for coadded wavelength At the moment, I assume a
    # linear dispersion for the data.  This is approximately correct
    # for the medium-resolution gratings, but definitely not correct
    # for G140L.  The G130M+G160M data get a dual-dispersion solution
    # broken at 1450A.

    goodwave = (fluxin != 0) & (wavein > 0.0)
    if chan == 3:
        # another case where the array should be float (double) type
        wave = disp[0]*np.array(range(33099), dtype=np.float64) + wavein[goodwave].min()
        good = np.where(wave <= 1460.0)
        wave = wave[good]
        wave = np.concatenate((wave,
                               np.array(wave.max() +
                                        disp[1] +
                                        disp[1]*np.array(range(27800),
                                                         dtype=np.float64))),
                              axis= 1)
        good = (wave >= wavein[goodwave].min()) & (wave <= wavein[goodwave].max())
        wave = wave[good]
    else:
        minwave = np.round(wavein[goodwave].min())
        maxwave = np.round(wavein[goodwave].max())
        wave = minwave + np.array(range(int((maxwave - minwave)/disp[chanind] + 1.0)),
                                  dtype=np.float64)*disp[chanind]

    # Interpolate the individual exposures onto the coadded wavelength
    # vector
        
    fluxint = np.zeros((len(wave), nfiles))
    errint = np.zeros((len(wave), nfiles))
    exptint = np.zeros((len(wave), nfiles))
    
    for i in range(nfiles):
        good = np.where(wavein[:, i] > wave.min())
        fluxtmp = np.interp(wave, wavein[good, i][0], fluxin[good, i][0], left = 0.0, right = 0.0)
        fluxint[:, i] = fluxtmp[:]
        errtmp = np.interp(wave, wavein[good, i][0], errin[good, i][0], left = 0.0, right = 0.0)
        errint[:, i] = errtmp[:]
        exptmp = np.interp(wave, wavein[good, i][0], exptin[good, i][0], left = 0.0, right = 0.0)
        exptint[:, i] = exptmp[:]

    # Perform the coaddition:
    flux = np.zeros(len(wave))
    err = np.zeros(len(wave))
    exptime = np.zeros(len(wave))

    for i in range(len(wave)):
        good = np.where(fluxint[i, :] != 0.0)
        good = good[0]
        if len(good) > 0:
            exptime[i] = sum(exptint[i, good])

            if len(good) == 1:
                flux[i] = fluxint[i, good][0]
                err[i] = errint[i, good][0]
            else:
                localsn = np.median(fluxint[i, good]/errint[i, good]) # local pixel S/N

                if method == -1:
                    flux[i] = sum(fluxint[i, good])/float(len(good))
                    err[i] = np.std(fluxint[i, good])/np.sqrt(float(len(good)))

                if method == 1:
                    flux[i] = sum(fluxint[i, good]*exptint[i, good])/sum(exptint[i, good])
                    err[i] = np.sqrt(sum((exptint[i, good]*errint[i, good])**2))/sum(exptint[i, good])

                if method == 2:
                    wts = (1.0/errint[i, good]**2)
                    flux[i] = sum(fluxint[i, good]*wts)/sum(wts)
                    err[i] = np.sqrt(1.0/sum(wts))

                if method == 3:
                    wts = (fluxint[i, good]/errint[i, good])**2
                    flux[i] = sum(fluxint[i, good]*wts)/sum(wts)
                    err[i] = np.sqrt(1.0/sum(wts))

    # end coaddition

    if ubershift != None:
        wave = wave + ubershift

    # Make cos spectrum object
    readme = "Generated by python-coaddx1d at %s" % datetime.now()
    spectrum = COSSpectrum(readme, targname, fitshdr, files, grating, side,
                           cenwave, fppos, dateobs, exptin, wavein, fluxin,
                           errin, wave, flux, err, exptime)

    # Save it one way or another...
    if savefile != None and savefile != False:
        if type(savefile) == str:
            fname = savefile
        elif type(savefile) == bool:
            fname = targname.lower() + "-" + channame[chan - 1].lower()
        else:
            raise Exception("keyword: savefile must be a string or boolean")
        
        if saveformat == "dat":
            f = open(fname + ".dat", 'wb')
            w = csv.writer(f, delimiter = ' ')
            try:
                for n, row in enumerate(wave):
                    w.writerow([row, flux[n], err[n], exptime[n]])
            except:
                f.close()
                raise Exception("Failure writing .dat file")
            f.close()
        elif saveformat == "csv":
            f = open(fname + ".csv", 'wb')
            w = csv.writer(f)
            try:
                for n, row in enumerate(wave):
                    w.writerow([row, flux[n], err[n], exptime[n]])
            except:
                f.close()
                raise Exception("Failure writing .csv file")
            f.close()
        elif saveformat == "pickle":
            f = open(fname + ".pkl", 'wb')
            pickle.dump(spectrum, f)
            f.close()
        else:
            # you're doing it wrong.
            if verbose:
                print("No valid format was specified. File not saved.")

    
    return spectrum
                

# Helper function that works as an IDL "clause"
def _herczeg(refx, refy, compx, compy, verbose = False):
    if abs(len(compx) - len(refx)) > 2:
        if verbose:
            if len(compx) < len(refx):
                # Pass all the data needed for the coaddx1d.pro print here...
                print("Comparison data does not cover the entire reference cross-correlation range")
                print("compx < refx")
            if len(refx) < len(compx):
                print("Reference data does not cover entire reference cross-correlation range.")
                print("refx < compx")

        refxrange = np.array([min(refx), max(refx)])
        compxrange = np.array([min(compx), max(compx)])
        userange = np.array([refxrange[0] if refxrange[0] > compxrange[0] else compxrange[0],
                             refxrange[1] if refxrange[1] > compxrange[1] else compxrange[1]])

        compgood = (compx >= userange[0]) & (compx <= userange[1])
        refgood = (refx >= userange[0]) & (refx <= userange[1])
        
        # compgood = np.concatenate((np.where(compx >= userange[0]),
        #                            np.where(compx <= userange[1])), axis = 1)
        # refgood = np.concatenate((np.where(refx >= userange[0]),
        #                           np.where(refx <= userange[1])), axis = 1)

        compx = compx[compgood]
        compy = compy[compgood]
        refx = refx[refgood]
        refy = refy[refgood]
        # RECURSION! 
        return _herczeg(refx, refy, compx, compy, verbose)
    else:
        return (refx, refy, compx, compy)
