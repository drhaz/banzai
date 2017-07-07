import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import dateutil.parser





def plotlongtermtrend (inputfile, instrument=None, filter=None):
    data = np.genfromtxt(inputfile, unpack=True, dtype=None,\
                  converters={ 1: lambda x: dateutil.parser.parse(x)}, names = ['name','dateobs', 'site', 'telescope', 'camera','filter','airmass','zp'])


    selection = np.ones(len (data['name']), dtype=bool)

    if filter is not None:
        selection = selection & (data['filter'] == filter )
    if instrument is not None:
        selection = selection & (data['camera'] == instrument)

    zpselect = data['zp'][selection]
    dateselect = data['dateobs'][selection]
    airmasselect = data['airmass'][selection]

    meanzp = np.nanmedian(zpselect)
    zp_air = zpselect + 0.2 * airmasselect - 0.2

    plt.plot (dateselect, zpselect, ".", c="grey", label="no airmass correction")
    plt.plot (dateselect, zp_air, ".", c="blue", label="with airmass correction" )
    plt.legend()
    plt.ylim([meanzp-0.5,meanzp+0.5])
    plt.gcf().autofmt_xdate()
    plt.xlabel ("DATE-OBS")
    plt.ylabel ("Photometric Zeropint %s" % (filter))
    plt.title ("Long term throughput  %s" % (instrument))
    plt.savefig ("photzptrend-%s-%s.png" % (instrument, filter))

    plt.figure()

    plt.plot (airmasselect, zpselect, ".", c="grey")
    plt.plot (airmasselect, zp_air, ".", c="blue")

    plt.xlabel ("Airmass")
    plt.ylabel ("Photomertic Zeropint %s" % (filter))
    plt.ylim([meanzp-0.5,meanzp+0.5])

    plt.savefig ("airmasstrend-%s-%s.png"  % (instrument, filter))


plotlongtermtrend ("photzp.db")