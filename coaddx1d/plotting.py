"""
Plotting functions
"""


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

    @type spectrum: L{COSSpectrum} object
    @var spectrum: a COSSpectrum object created by coaddx1d.
    @type err: bool
    @var err: boolean switch for plotting the flux error (default True)
    @type scale: float
    @var scale: scale*mean(flux) determines the flux axis height (default 2.0)
                This is a quick and dirty way of adjusting the flux axis upper
                limit.

    @return fig: a matplotlib figure instance of the figure
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

    if spectrum.grating[0] == "G130M" or \
       spectrum.grating[0] == "G160M" or \
       spectrum.grating[0] == "G140L":
        panes = 2
    else:
        # G130M + G160M
        panes = 4
    

    panesize = (spectrum.wave.max() - spectrum.wave.min() + 2.0*buff)/2.0
    for i in range(panes):
        ax = mp.subplot(2, 1, (i + 1))
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

    @type spectrum: L{COSSpectrum object}
    @param spectrum: a COSSpectrum object created by coaddx1d

    @return fig: a matplotlib figure instance
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
