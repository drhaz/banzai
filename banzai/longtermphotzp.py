
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
import dateutil.parser
import sys
import calendar


airmasscorrection = {'gp':0.2, 'rp': 0.12, 'ip': 0.08, 'zp': 0.05,}

def readDataFile (inputfile):
    return np.genfromtxt(inputfile, unpack=True, dtype=None,skip_footer=5, \
                    converters={ 1: lambda x: dateutil.parser.parse(x)}, names = ['name','dateobs', 'site', 'dome', 'telescope', 'camera','filter','airmass','zp'])


def findUpperEnvelope (dateobs, datum, range=1, ymax = 24.2):
    """
    Find the upper envelope of a photZP time line

    Idea:
    For a set of day(s), find the highest n zeropoint measurements within an error range.
    Omit the brightest one for outlier rejection. Average the remaining data points of the day. Accept as value

    After this, filter value. Currently, I use a cheap-o-Kalman fake with a fixed Kalman gain. Not smart, but works

    Reject some bright and too faint outliers, i.e., if Kalman correction is too large, reject as an obvious large issue.


    :param dateobs:
    :param datum:
    :param range:
    :return:
    """

    stderror = 0.02

    alldata = zip (dateobs, datum)
    sorted_points = sorted (alldata)
    x = np.asarray([point[0] for point in sorted_points])
    y = np.asarray([point[1] for point in sorted_points])

    day_x = []
    day_y = []

    # define  day / night boundary for the site.
    startdate =datetime.datetime(year=x[0].year, month = x[0].month, day=x[0].day, hour=20)
    enddate = x[len(x)-1]
    while startdate < enddate:
        todayzps = y [ (x > startdate ) & (x < startdate + datetime.timedelta(days=3)) & (y<ymax) & (y is not np.nan)]

        if len (todayzps) > 3: # require a minium amount of data for a night
            todayzps = np.sort (todayzps)[1:]
            max = np.nanmax(todayzps)
            upperEnv = np.nanmean(todayzps [ todayzps > (max - stderror)])
            if upperEnv is not np.nan:
                day_x.append (startdate)
                day_y.append (upperEnv)

        startdate = startdate + datetime.timedelta(days=1)


    # filter the daily zero point variation. Work in progress.

    newday_y = [day_y[0]]

    for y in day_y:

        last = newday_y[len (newday_y)-1]
        correction = y - last

        if (np.abs(correction) < 0.1):
            upd = last + 1* correction
        elif np.abs (correction) < 0.4:
            upd = last + 0.2 * correction
        else:
            upd = last
        newday_y.append (upd)

    return np.asarray(day_x), np.asarray(newday_y[1:])


def trendcorrectthroughput (datadate, datazp, modeldate, modelzp):

    """ Detrend input data based in a trend model

    """

    modelgmt = np.zeros((len (modeldate)))
    for ii in range(0,len(modeldate)):
        modelgmt[ii] = calendar.timegm(modeldate[ii].timetuple())

    corrected = np.zeros(len(datazp))
    for ii in range (0, len(corrected)):
        interpolated = np.interp(calendar.timegm(datadate[ii].timetuple()), modelgmt, modelzp)
        corrected[ii] = datazp[ii] - interpolated

    return corrected

def plotallmirrormodels (basedirectory="/home/dharbeck/lcozpplots"):
    import glob
    modellist = glob.glob ("%s/mirrormodel*gp.dat" % (basedirectory))

    for model in modellist:
        data = np.genfromtxt (model, dtype=None, names=("date", "time", "zp") )
        datestring = np.core.defchararray.add(data['date'] , "T" )
        datestring = np.core.defchararray.add(datestring, data['time'] )
        date = [np.datetime64(x) for x in datestring]
        date = np.asarray(date).astype(datetime.datetime)

        #data['zp'] = np.power(10, data['zp']/2.5)
        #data['zp'] = data['zp'] / np.min(data['zp'])


        plt.gcf().autofmt_xdate()
        plt.plot (date, data['zp'], label=model[-11:-7])


        # ax = fig.gca()
        # ax.set_xticks(numpy.arange(0, 1, 0.1))
        # ax.set_yticks(numpy.arange(0, 1., 0.1))
    plt.legend()

    plt.grid(True, which='both')
    plt.savefig ("%s/allmodels.png" % basedirectory)
    plt.close()

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

    ymax = 24.2  # good starting point for 2m:spectral cameras
    if (instrument is not None):
        if instrument.startswith("fl"): # 1m sinistro
            ymax = 23.75
        if instrument.startswith("kb"): # 0.4m sbigs
            ymax = 22

    meanzp = np.nanmedian(zpselect)
    zp_air = zpselect + airmasscorrection[filter] * airmasselect - airmasscorrection[filter]

    # find the overall trend of zeropoint variations.
    _x, _y = findUpperEnvelope(dateselect, zp_air, ymax= ymax)
    detrended = trendcorrectthroughput(dateselect, zp_air, _x, _y)

    plt.figure()
    #plt.plot (dateselect, zpselect, ".", c="grey", label="no airmass correction")
    plt.plot (dateselect,  zp_air, ".", c="blue", label="with airmass correction" )
    plt.plot (_x, _y, "-", c='red', label='upper envelope')

    plt.plot (dateselect, detrended+ymax - 1.5,  ".", c="cyan", label="detrended + 23mag")
    plt.legend()
    plt.xlim ([datetime.datetime(2016,01,01),datetime.datetime (2017,9,01)])
    plt.ylim([ymax - 2,ymax])
    plt.gcf().autofmt_xdate()
    plt.xlabel ("DATE-OBS")
    plt.ylabel ("Photometric Zeropint %s" % (filter))
    plt.title ("Long term throughput  %s in %s" % (instrument,filter))

    if (instrument in ("kb97","kb98")):
        plt.axvline(x=datetime.datetime(2017, 06,30), color='k', linestyle='--')

    outfigname = "%s/photzptrend-%s-%s.png" % (basedirectory,instrument, filter)
    plt.savefig (outfigname)
    plt.close()

    outdetrendfname = "%s/photdetrend-%s-%s.dat" % (basedirectory,instrument, filter)

    np.savetxt(outdetrendfname, np.c_[dateselect, zp_air, detrended], header="DATE-OBS photzp photzp_detrended",
               fmt="%s %f %f")

    outmodelfname = "%s/mirrormodel-%s-%s.dat" % (basedirectory,instrument, filter)

    np.savetxt(outmodelfname, np.c_[_x, _y], header="DATE-OBS zp envelope",
           fmt="%s %f ")

    plt.close()
    plt.figure()

    plt.plot (airmasselect, zpselect, ".", c="grey")
    plt.plot (airmasselect, zp_air, ".", c="blue")

    plt.xlabel ("Airmass")
    plt.ylabel ("Photomertic Zeropint %s" % (filter))
    plt.ylim([meanzp-0.5,meanzp+0.5])

    plt.savefig ("%s/airmasstrend-%s-%s.png"  % (basedirectory, instrument, filter))
    plt.close()


if __name__ == '__main__':

    plotallmirrormodels()
    sys.exit(0)

    plotlongtermtrend ("ogg", filter='gp', instrument='fs02')
    plotlongtermtrend ("ogg", filter='rp', instrument='fs02')
    #
    plotlongtermtrend ("coj", filter='gp', instrument='fs01')
    plotlongtermtrend ("coj", filter='rp', instrument='fs01')

    plotlongtermtrend ("elp", filter='gp', instrument='fl05')
    plotlongtermtrend ("elp", filter='rp', instrument='fl05')

  #  sys.exit(1)
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

