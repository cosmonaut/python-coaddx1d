"""
Plotting functions

Author: Nico Nell (nicholas.nell@colorado.edu)
"""


import numpy as np
from scipy.signal import medfilt

try:
    import matplotlib.pyplot as mp
    from matplotlib.ticker import ScalarFormatter
    from matplotlib.ticker import MultipleLocator
except:
    print("WARNING: matplotlib was not imported. Plotting functions will not work")


def plotflux(spectrum, err = True, scale = 2.0):
    """
    Create a flux vs. wavelength plot.

    Require Arguments:

        *spectrum*: [ COSSpectrum ]
        a COSSpectrum object created by coaddx1d.
        
        *err*: [ bool ]
        boolean switch for plotting the flux error (default True)
        
        *scale*: [ float ]
        scale*mean(flux) determines the flux axis height (default
        2.0).  This is a quick and convenient way of adjusting the
        flux axis upper limit.

    Returns:
    
        *fig*: [ matplotlib.pyplot.figure ]
        a matplotlib figure instance of the figure
    """

    fig = mp.figure()
    yfmt = ScalarFormatter(useMathText = True)
    
    xlabel = r'Wavelength ($\AA$)'
    ylabel = r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)'

    buff = 10.0

    xmin = spectrum.wave.min() - buff
    xmax = spectrum.wave.max() + buff
    ymin = 0.0
    ymax = spectrum.flux.mean()*scale

    mp.ylim((ymin, ymax))

    # Parse grating array...
    
    if len(np.unique(spectrum.grating)) > 1:
        # G130M" + "G160M
        panes = 4
    else:
        # Just one grating 
        panes = 2

    panesize = (spectrum.wave.max() - spectrum.wave.min() + 2.0*buff)/panes
    
    for i in range(panes):
        ax = mp.subplot(panes, 1, (i + 1))
        ax.yaxis.set_major_formatter(yfmt)
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(5))

        xmin = (spectrum.wave.min() - buff) + panesize*float(i)
        xmax = (spectrum.wave.min() - buff) + panesize*float(i + 1.0)
        
        mp.xlim((xmin, xmax))
        mp.ylim((ymin, ymax))
        mp.xlabel(xlabel)
        mp.ylabel(ylabel)

        region = (spectrum.wave >= xmin) & (spectrum.wave <= xmax)
        
        mp.plot(spectrum.wave[region],
                medfilt(spectrum.flux[region], 7.0),
                color = 'black',
                linewidth = 0.5)
        
        if err:
            mp.plot(spectrum.wave[region],
                    medfilt(spectrum.err[region], 7.0),
                    color = 'red',
                    linewidth = 0.5,
                    linestyle = 'dotted')

    
    # if show:
    #     mp.show()
    # else:
    #     #mp.savefig(name)

    return fig


def plotexptime(spectrum, show = False):
    """
    Create an exposure time vs. wavelength plot

    Require Arguments:

        *spectrum*: [ COSSpectrum ]
        a COSSpectrum object created by coaddx1d.

    Returns:
    
        *fig*: [ matplotlib.pyplot.figure ]
        a matplotlib figure instance of the figure
    """

    # Just a basic implementation for now...
    fig = mp.figure()

    xlabel = r'Wavelength ($\AA$)'
    ylabel = r'Exposure Time (s)'

    mp.xlabel(xlabel)
    mp.ylabel(ylabel)

    # We should divide this for G130M + G160M?
    mp.plot(spectrum.wave,
            medfilt(spectrum.exptime, 7.0),
            color = 'black',
            linewidth = 0.5)

    return fig
