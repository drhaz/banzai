
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
import dateutil.parser


airmasscorrection = {'gp':0.2, 'rp': 0.12, 'ip': 0.08, 'zp': 0.05,}

def readDataFile (inputfile):
    return np.genfromtxt(inputfile, unpack=True, dtype=None,skip_footer=5, \
                    converters={ 1: lambda x: dateutil.parser.parse(x)}, names = ['name','dateobs', 'site', 'dome', 'telescope', 'camera','filter','airmass','zp'])


def plotlongtermtrend (site, enclosure=None, telescope=None, instrument=None, filter=None, basedirectory="/home/dharbeck/lcozpplots"):

    inputfile = "%s/%s-%s.db" % (basedirectory, site, instrument)

    mirrorfilename = "%s/mirror_%s.db" % (basedirectory,instrument)
    data = readDataFile(inputfile)

    selection = np.ones(len (data['name']), dtype=bool)

    if filter is not None:
        selection = selection & (data['filter'] == filter )
    if instrument is not None:
        selection = selection & (data['camera'] == instrument)

    selection = selection & np.logical_not (  np.isnan(data['zp']))

    zpselect = data['zp'][selection]
    dateselect = data['dateobs'][selection]
    airmasselect = data['airmass'][selection]



    ymax = 24.2
    if (instrument is not None):
        if instrument.startswith("fl"):
            ymax = 23.75
        if instrument.startswith("kb"):
            ymax = 22

    meanzp = np.nanmedian(zpselect)
    zp_air = zpselect + airmasscorrection[filter] * airmasselect - airmasscorrection[filter]


    plt.figure()
    #plt.plot (dateselect, zpselect, ".", c="grey", label="no airmass correction")
    plt.plot (dateselect, zp_air, ".", c="blue", label="with airmass correction" )
    plt.legend()
    plt.xlim ([datetime.datetime(2016,01,01),datetime.datetime (2017,8,01)])
    plt.ylim([ymax - 2,ymax])
    plt.gcf().autofmt_xdate()
    plt.xlabel ("DATE-OBS")
    plt.ylabel ("Photometric Zeropint %s" % (filter))
    plt.title ("Long term throughput  %s in %s" % (instrument,filter))

    if (instrument in ("kb97","kb98")):
        plt.axvline(x=datetime.datetime(2017, 06,30), color='k', linestyle='--')

    plt.savefig ("%s/photzptrend-%s-%s.png" % (basedirectory,instrument, filter))
    plt.close()
    plt.figure()

    plt.plot (airmasselect, zpselect, ".", c="grey")
    plt.plot (airmasselect, zp_air, ".", c="blue")

    plt.xlabel ("Airmass")
    plt.ylabel ("Photomertic Zeropint %s" % (filter))
    plt.ylim([meanzp-0.5,meanzp+0.5])

    plt.savefig ("%s/airmasstrend-%s-%s.png"  % (basedirectory, instrument, filter))



if __name__ == '__main__':

    plotlongtermtrend ("ogg", filter='gp', instrument='fs02')
    plotlongtermtrend ("ogg", filter='rp', instrument='fs02')
    #
    plotlongtermtrend ("coj", filter='gp', instrument='fs01')
    plotlongtermtrend ("coj", filter='rp', instrument='fs01')

    plotlongtermtrend ("elp", filter='gp', instrument='fl05')
    plotlongtermtrend ("elp", filter='rp', instrument='fl05')
#
#
    plotlongtermtrend ("lsc", filter='gp', instrument='fl03')
    plotlongtermtrend ("lsc", filter='rp', instrument='fl03')
#
    plotlongtermtrend ("lsc", filter='gp', instrument='fl04')
    plotlongtermtrend ("lsc", filter='rp', instrument='fl04')
#
    plotlongtermtrend ("cpt", filter='gp', instrument='fl06')
    plotlongtermtrend ("cpt", filter='rp', instrument='fl06')
#
#
    plotlongtermtrend ("lsc", filter='gp', instrument='fl15')
    plotlongtermtrend ("lsc", filter='rp', instrument='fl15')
#
    plotlongtermtrend ("cpt", filter='gp', instrument='fl14')
    plotlongtermtrend ("cpt", filter='rp', instrument='fl14')
#
#
    plotlongtermtrend ("cpt", filter='gp', instrument='fl16')
    plotlongtermtrend ("cpt", filter='rp', instrument='fl16')

    plotlongtermtrend ("coj", filter='gp', instrument='kb97')
    plotlongtermtrend ("coj", filter='rp', instrument='kb97')
    #
    plotlongtermtrend ("tfn", filter='gp', instrument='kb29')
    plotlongtermtrend ("tfn", filter='rp', instrument='kb29')
    #
    plotlongtermtrend ("ogg", filter='gp', instrument='kb27')
    plotlongtermtrend ("ogg", filter='rp', instrument='kb27')
    #
    plotlongtermtrend ("ogg", filter='gp', instrument='kb82')
    plotlongtermtrend ("ogg", filter='rp', instrument='kb82')
    #
    #
    #
    plotlongtermtrend ("coj", filter='gp', instrument='kb98')
    plotlongtermtrend ("coj", filter='rp', instrument='kb98')
