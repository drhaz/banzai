import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
import dateutil.parser
import sys
import calendar
from StringIO import StringIO

airmasscorrection = {'gp': 0.17, 'rp': 0.09, 'ip': 0.06, 'zp': 0.05, }

colorterms = {}


def readDataFile(inputfile):
    file = open(inputfile)
    contents = file.read()
    contents = contents.replace('UNKNOWN', 'nan')
    file.close()
    data = np.genfromtxt(StringIO(contents), unpack=True, dtype=None,
                         skip_footer=5, \
                         converters={1: lambda x: dateutil.parser.parse(x)},
                         names=['name', 'dateobs', 'site', 'dome',
                                'telescope', 'camera', 'filter', 'airmass',
                                'zp', 'colorterm', 'zpsig'])

    return data


def findUpperEnvelope(dateobs, datum, range=1, ymax=24.2):
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

    alldata = zip(dateobs, datum)
    sorted_points = sorted(alldata)
    x = np.asarray([point[0] for point in sorted_points])
    y = np.asarray([point[1] for point in sorted_points])

    day_x = []
    day_y = []

    # TODO: define  day / night boundary for the site.
    startdate = datetime.datetime(year=x[0].year, month=x[0].month,
                                  day=x[0].day, hour=12)
    enddate = x[len(x) - 1]
    while startdate < enddate:
        todayzps = y[
            (x > startdate) & (x < startdate + datetime.timedelta(days=10)) & (
                y < ymax) & (y is not np.nan)]

        if len(todayzps) > 5:  # require a minium amount of data for a night
            todayzps = np.sort(todayzps)[3:]
            max = np.nanmax(todayzps)
            upperEnv = np.nanmean(todayzps[todayzps > (max - stderror)])
            if upperEnv is not np.nan:
                day_x.append(startdate)
                day_y.append(upperEnv)

        startdate = startdate + datetime.timedelta(days=1)

    # filter the daily zero point variation. Work in progress.


    if len(day_y) > 0:

        newday_y = [day_y[0]]

        for y in day_y:

            last = newday_y[len(newday_y) - 1]
            correction = y - last

            if (np.abs(correction) < 0.1):
                upd = last + 1 * correction
            elif np.abs(correction) < 0.4:
                upd = last + 0.2 * correction
            else:
                upd = last + 0.1 * correction
            newday_y.append(upd)
    else:
        newday_y = [0, 0]

    return np.asarray(day_x), np.asarray(newday_y[1:])


def trendcorrectthroughput(datadate, datazp, modeldate, modelzp):
    """ Detrend input data based in a trend model

    """

    modelgmt = np.zeros((len(modeldate)))
    for ii in range(0, len(modeldate)):
        modelgmt[ii] = calendar.timegm(modeldate[ii].timetuple())

    corrected = np.zeros(len(datazp))
    for ii in range(0, len(corrected)):
        interpolated = np.interp(calendar.timegm(datadate[ii].timetuple()),
                                 modelgmt, modelzp)
        corrected[ii] = datazp[ii] - interpolated

    # estimate if photometric


    photomerticthres = 0.25
    day_x = []
    day_y = []
    alldata = zip(datadate, corrected)
    sorted_points = sorted(alldata)
    x = np.asarray([point[0] for point in sorted_points])
    y = np.asarray([point[1] for point in sorted_points])
    startdate = datetime.datetime(year=2016, month=4, day=1, hour=16)
    enddate = x[len(x) - 1]

    while startdate < enddate:
        todayzps = y[
            (x > startdate) & (x < startdate + datetime.timedelta(days=1))]

        photometric = -1

        if len(todayzps) > 0:  # require a minium amount of data for a night

            if np.min(todayzps > -0.15):
                photometric = 1
            else:
                photometric = 0

        day_x.append(startdate)
        day_y.append(photometric)

        startdate = startdate + datetime.timedelta(days=1)

    day_x = np.asarray(day_x)
    day_y = np.asarray(day_y)
    unclassified = len(day_y[day_y < 0])
    photometric = len(day_y[day_y > 0])
    nonphot = len(day_y[day_y == 0])
    #

    # print ("out of %d days\nphotometric\t%d\nnon-photometric\t%d\nunknown\t%d" %
    #        (unclassified+photometric+nonphot, photometric, nonphot, unclassified))

    return corrected, day_x, day_y


def plotallmirrormodels(basedirectory="/home/dharbeck/lcozpplots"):
    import glob
    modellist = glob.glob("%s/mirrormodel-[fl|fs]*rp.dat" % (basedirectory))

    for model in modellist:
        print model
        data = np.genfromtxt(model, dtype=None, names=("date", "time", "zp"))
        datestring = np.core.defchararray.add(data['date'], "T")
        datestring = np.core.defchararray.add(datestring, data['time'])
        date = [np.datetime64(x) for x in datestring]
        date = np.asarray(date).astype(datetime.datetime)

        # data['zp'] = np.power(10, data['zp']/2.5)
        # data['zp'] = data['zp'] / np.min(data['zp'])

        plt.gcf().autofmt_xdate()
        plt.plot(date, data['zp'], label=model[-11:-7])

        # ax = fig.gca()
        # ax.set_xticks(numpy.arange(0, 1, 0.1))
        # ax.set_yticks(numpy.arange(0, 1., 0.1))
    # plt.legend()
    plt.ylabel("phot zeropoint rp")
    plt.xlim([datetime.datetime(2016, 01, 01), datetime.datetime(2017, 11, 01)])
    plt.grid(True, which='both')
    plt.savefig("%s/allmodels.png" % basedirectory)
    plt.close()


def plotlongtermtrend(site, enclosure=None, telescope=None, instrument=None,
                      filter=None, basedirectory="/home/dharbeck/lcozpplots"):
    print site, telescope, instrument
    inputfile = "%s/%s-%s.db" % (basedirectory, site, instrument)
    mirrorfilename = "%s/mirror_%s.db" % (basedirectory, instrument)
    data = readDataFile(inputfile)



    # down-select data by viability and camera / filer combination
    selection = np.ones(len(data['name']), dtype=bool)

    if filter is not None:
        selection = selection & (data['filter'] == filter)
    if instrument is not None:
        selection = selection & (data['camera'] == instrument)
    selection = selection & np.logical_not(np.isnan(data['zp']))
    selection = selection & np.logical_not(np.isnan(data['airmass']))

    if (len(selection) == 0):
        print (
            "Zero viable elements left for %s %s. Ignoring" % (
            site, instrument))
        return

    zpselect = data['zp'][selection]
    dateselect = data['dateobs'][selection]
    airmasselect = data['airmass'][selection]
    try:
        zpsigselect = data['zpsig'][selection]
    except:
        zpsigselect = zpsigselect * 0
    ymax = 24.2  # good starting point for 2m:spectral cameras
    photzpmaxnoise = 0.2
    if (instrument is not None):
        if instrument.startswith("fl"):  # 1m sinistro
            ymax = 23.95
        if instrument.startswith("kb"):  # 0.4m sbigs
            ymax = 22
            photzpmaxnoise = 0.5

    # Calculate air-mass corrected photometric zeropoint
    zp_air = zpselect + airmasscorrection[filter] * airmasselect - \
             airmasscorrection[filter]

    # find the overall trend of zeropoint variations.
    detrend = photdate = photflat = None
    try:
        _x, _y = findUpperEnvelope(dateselect[zpsigselect<photzpmaxnoise], zp_air[zpsigselect<photzpmaxnoise],
                                   ymax=ymax)

        outmodelfname = "%s/mirrormodel-%s-%s.dat" % (
            basedirectory, instrument, filter)
        np.savetxt(outmodelfname, np.c_[_x, _y], header="DATE-OBS zp envelope",
                   fmt="%s %f ")

        detrended = _y
        # detrended, photdate, photflag = trendcorrectthroughput(dateselect[zpsigselect<photzpmaxnoise],
        #                                                           zp_air[zpsigselect<photzpmaxnoise], _x, _y)
    except Exception as e:
        detrend = None
        _x = None
        print ("Failure while detrending!", e)
    plt.figure()
    # plt.plot (dateselect, zpselect, ".", c="grey", label="no airmass correction")

    plt.plot(dateselect[zpsigselect >= photzpmaxnoise], zp_air[zpsigselect >= photzpmaxnoise],
             '.',
             markersize=1,
             c="grey",
             label="large error")
    plt.plot(dateselect[zpsigselect < photzpmaxnoise], zp_air[zpsigselect < photzpmaxnoise], 'o',
            
             markersize=2,
             c="blue",
             label="small error")

    if _x is not None:
        plt.plot(_x, _y, "-", c='red', label='upper envelope')
        # plt.plot(dateselect, detrended + ymax - 1.5, ".", c="cyan",
        #          label="detrended + 23mag")
        #
        # plt.plot (photdate, photflag/10. + ymax - 1.5, ".", c='green', label = "photometric flag")

    plt.legend()
    plt.xlim([datetime.datetime(2016, 01, 01), datetime.datetime(2017, 11, 01)])
    plt.ylim([ymax - 2, ymax])
    plt.gcf().autofmt_xdate()
    plt.xlabel("DATE-OBS")
    plt.ylabel("Photometric Zeropoint %s" % (filter))
    plt.title("Long term throughput  %s in %s" % (instrument, filter))

    if (instrument in ("kb97", "kb98")):
        plt.axvline(x=datetime.datetime(2017, 6, 30), color='k', linestyle='--')
    if (instrument in ("fl03", "xxx")):
        plt.axvline(x=datetime.datetime(2017, 8, 31), color='k', linestyle='--')

    outfigname = "%s/photzptrend-%s-%s.png" % (
        basedirectory, instrument, filter)
    plt.savefig(outfigname, dpi=600)
    plt.close()

    outdetrendfname = "%s/photdetrend-%s-%s.dat" % (
        basedirectory, instrument, filter)

    if detrend is not None:
        np.savetxt(outdetrendfname, np.c_[dateselect, zp_air, detrended],
                   header="DATE-OBS photzp photzp_detrended",
                   fmt="%s %f %f")

    plt.close()

    plt.figure()

    plt.hist(zpsigselect, 50, range=[0, 1], normed=True)
    outerrorhistfname = "%s/errorhist-%s-%s.png" % (
        basedirectory, instrument, filter)
    plt.savefig(outerrorhistfname)
    plt.close()

    plt.figure()

    plt.plot(airmasselect, zpselect, ".", c="grey")
    plt.plot(airmasselect, zp_air, ".", c="blue")

    plt.xlabel("Airmass")
    plt.ylabel("Photomertic Zeropoint %s" % (filter))
    meanzp = np.nanmedian(zpselect)
    plt.ylim([meanzp - 0.5, meanzp + 0.5])

    plt.savefig(
        "%s/airmasstrend-%s-%s.png" % (basedirectory, instrument, filter))
    plt.close()



    ### Color terms
    plt.figure()
    selection = selection & np.logical_not(np.isnan(data['colorterm']))
    selection = selection & (np.abs(data['colorterm']) < 0.3)
    selection_lonoise = selection & (data['zpsig'] < 0.2)
    selection_hinoise = selection & (data['zpsig'] >= 0.2)

    plt.plot(data['dateobs'][selection_hinoise], data['colorterm'][
        selection_hinoise], '.', markersize=2, c="grey",
             label="color term [ hi sigma] %s " % (filter))


    colortermselect = data['colorterm'][selection_lonoise]
    dateselect = data['dateobs'][selection_lonoise]
    meancolorterm = np.median(colortermselect)
    plt.plot(dateselect, colortermselect, 'o', markersize=2, c="blue",
             label="color term [low sigma] %s " % (filter))
    plt.axhline(y=meancolorterm, color='r', linestyle='-')
    print ("Color term filter %s : % 5.3f" % (filter, meancolorterm))

    if filter not in colorterms:
        colorterms[filter] = {}
    colorterms[filter][instrument] = meancolorterm

    plt.xlim([datetime.datetime(2016, 01, 01), datetime.datetime(2017, 11, 01)])
    plt.ylim([-0.2, 0.2])

    plt.savefig(
        "%s/colortermtrend-%s-%s.png" % (basedirectory, instrument, filter))
    plt.close()


def plotallcolorterms():
    for filter in colorterms:

        for camera in colorterms[filter]:
            print (
                "% 3s % 4s % 5.2f" % (
                filter, camera, colorterms[filter][camera]))

        cameras = colorterms[filter]
        terms = []
        for camera in cameras:
            terms.append(colorterms[filter][camera])
        plt.plot(cameras, terms, "o")
        plt.title(filter)
        plt.savefig("colorterms-%s.png" % (filter))
        plt.close()
    pass


import os
import re

if __name__ == '__main__':
    plt.style.use('ggplot')
    basedirectory = "/Users/harbeck/lcoplots"
    databases = [each for each in os.listdir(basedirectory) if
                 each.endswith('.db')]

    for db in databases:
        match = re.search('(\D\D\D)-(\D\D\d\d)\.db', db)
        if match:
            for filter in ( 'gp', 'rp'):
                site = match.group(1)
                camera = match.group(2)

                plotlongtermtrend(site, filter=filter, instrument=camera,
                                       basedirectory=basedirectory)

    plotallmirrormodels(basedirectory = basedirectory)

    #plotallcolorterms ()

    sys.exit(0)
