import matplotlib.pyplot as mp
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
#from matplotlib.ticker import LinearLocator
from scipy.signal import medfilt

def plotflux(spectrum, err = True, show = False, name = None):
    """
    Create a flux vs. wavelength plot.

    Inputs:
    spectrum: a COSSpectrum object created by coaddx1d.
    err: boolean switch for plotting the flux error (default True)

    Output:
    fig: a matplotlib figure instance of the figure
    """

    mp.clf()

    fig = mp.figure()
    yfmt = ScalarFormatter(useMathText = True)
    
    xlabel = r'Wavelength ($\AA$)'
    ylabel = r'Flux (erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)'

    buff = 10.0

    xmin = spectrum.wave.min() - buff
    xmax = spectrum.wave.max() + buff
    ymin = 0.0
    ymax = spectrum.flux.mean()*2.0

    mp.ylim((ymin, ymax))

    if spectrum.grating[0] == "G130M":
        half = (spectrum.wave.max() - spectrum.wave.min() + 2.0*buff)/2.0
        for i in range(2):
            ax = mp.subplot(2, 1, (i + 1)%2)
            ax.yaxis.set_major_formatter(yfmt)
            ax.xaxis.set_major_locator(MultipleLocator(20))
            ax.xaxis.set_minor_locator(MultipleLocator(5))

            xmin = (spectrum.wave.min() - buff) + half*float(i)
            xmax = (spectrum.wave.min() - buff) + half*float(i + 1.0)

            mp.xlim((xmin, xmax))
            mp.ylim((ymin, ymax))
            mp.xlabel(xlabel)
            mp.ylabel(ylabel)

            mp.plot(spectrum.wave,
                    medfilt(spectrum.flux, 7.0),
                    color = 'black',
                    linewidth = 0.5)

            if err:
                mp.plot(spectrum.wave,
                        medfilt(spectrum.err, 7.0),
                        color = 'red',
                        linewidth = 0.5,
                        linestyle = 'dotted')

    
    # if show:
    #     mp.show()
    # else:
    #     #mp.savefig(name)

    #mp.clf()
    return fig


def plotexptime(spectrum, show = False):
    """
    Create an exposure time vs. wavelength plot

    Inputs:
    spectrum: a COSSpectrum object created by coaddx1d

    Output:
    fig: a matplotlib figure instance
    """

    # Just a basic implementation for now...
    mp.clf()
    fig = mp.figure()

    # We should divide this for G130M + G160M?
    mp.plot(spectrum.wave,
            spectrum.exptime,
            color = 'black',
            linewidth = 0.5)

    return fig
