from __future__ import absolute_import, division, print_function, unicode_literals
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
import sys
import calendar
import astropy.time as astt
import scipy.signal
import argparse
import logging
import glob
import os
from astropy.io import ascii



airmasscorrection = {'gp': 0.17, 'rp': 0.09, 'ip': 0.06, 'zp': 0.05, }

colorterms = {}
telescopedict = {
    'lsc': ['doma:1m0a', 'domb:1m0a', 'domc:1m0a', 'aqwa:0m4a', 'aqwb:0m4a'],
    'coj': ['clma:2m0a', 'doma:1m0a', 'domb:1m0a', 'clma:0m4a', 'clma:0m4b'],
    'ogg': ['clma:2m0a', 'clma:0m4a', 'clma:0m4b'],
    'elp': ['doma:1m0a', 'aqwa:0m4a'],
    'cpt': ['doma:1m0a', 'domb:1m0a', 'domc:1m0a'],
    'tfn': ['aqwa:0m4a', 'aqwa:0m4b'],
    'sqa': ['doma:0m8a']
}


def readDataFile(inputfile):
    with open(inputfile, 'r') as file:
        contents = file.read()
        file.close()

        contents = contents.replace('UNKNOWN', 'nan')
        data = ascii.read(contents, names=['name', 'dateobs', 'site', 'dome',
                                           'telescope', 'camera', 'filter', 'airmass',
                                           'zp', 'colorterm', 'zpsig'], )

        data['dateobs'] = astt.Time(data['dateobs'], scale='utc', format='isot').to_datetime()

        return data

    return None


def getCombineddataByTelescope(site, telescope, context, instrument=None):
    """
    Concatenate all zeropoint data for a site, and select by telescope and instrument.
    :param site:
    :param telescope: string slecting dome *& telescope: 'domb:1m0a'
    :param context:
    :param instrument:
    :return: concatenated data for a site / tel / isntrument selection.
    """

    inputfiles = glob.glob("%s/%s-%s.db" % (context.imagedbPrefix, site, '*' if instrument is None else instrument))
    alldata = None

    for inputfile in inputfiles:
        # print ('Reading in %s' % inputfile)
        data = readDataFile(inputfile)
        if data is None:
            continue
        if alldata is None:
            alldata = data
        else:
            try:
                alldata = np.append(alldata, data)
            except Exception as e:
                print("Failed to append data for file %s" % inputfile, e)

    if alldata is None:
        return None

    # Helpful diagnostic tool to see what data are pulled in.
    # domedict = np.unique (alldata['dome'])
    # for dome in domedict:
    #     teldict = np.unique(alldata [ alldata['dome'] == dome] ['telescope'])
    #     print (site, dome, teldict)

    dome, tel = telescope.split(':')
    selection = (alldata['dome'] == dome) & (alldata['telescope'] == tel)
    alldata = alldata[selection]
    return alldata


def plotlongtermtrend(select_site, select_telescope, select_filter, context, instrument=None):
    data = getCombineddataByTelescope(select_site, select_telescope, context, instrument)

    if data is None:
        return
    # down-select data by viability and camera / filer combination
    selection = np.ones(len(data['name']), dtype=bool)

    if select_filter is not None:
        selection = selection & (data['filter'] == select_filter)
    if instrument is not None:
        selection = selection & (data['camera'] == instrument)

    # weed out bad data
    selection = selection & np.logical_not(np.isnan(data['zp']))
    selection = selection & np.logical_not(np.isnan(data['airmass']))

    if len(selection) == 0:
        return

    zpselect = data['zp'][selection]
    dateselect = data['dateobs'][selection]
    airmasselect = data['airmass'][selection]
    cameraselect = data['camera'][selection]
    zpsigselect = data['zpsig'][selection]

    ymax = 25.5  # good starting point for 2m:spectral cameras
    photzpmaxnoise = 0.2
    if select_telescope is not None:

        if '0m4' in select_telescope:  # 0.4m sbigs
            ymax = 22.5
            photzpmaxnoise = 0.5

    # Calculate air-mass corrected photometric zeropoint
    zp_air = zpselect + airmasscorrection[select_filter] * airmasselect - airmasscorrection[select_filter]

    # find the overall trend of zeropoint variations, save to output file.
    _x, _y = findUpperEnvelope(dateselect[zpsigselect < photzpmaxnoise], zp_air[zpsigselect < photzpmaxnoise],
                               ymax=ymax)
    outmodelfname = "%s/mirrormodel-%s-%s-%s.dat" % (
        context.imagedbPrefix, select_site, select_telescope, select_filter)
    np.savetxt(outmodelfname, np.c_[_x, _y], header="DATE-OBS zp envelope",
               fmt="%s %f ")

    plt.figure()

    uniquecameras = np.unique(cameraselect)
    for uc in uniquecameras:
        # plot zeropoint with differnt markers per camera

        plt.plot(dateselect[(zpsigselect <= photzpmaxnoise) & (cameraselect == uc)],
                 zp_air[(zpsigselect <= photzpmaxnoise) & (cameraselect == uc)],
                 'o', markersize=2, label=uc)
        plt.plot(dateselect[zpsigselect > photzpmaxnoise], zp_air[zpsigselect > photzpmaxnoise], '.',
                 markersize=1, c="grey", )

    if _x is not None:
        plt.plot(_x, _y, "-", c='red', label='upper envelope')

    else:
        print("Mirror model failed to compute. not plotting !")

    plt.legend()
    plt.xlim([datetime.datetime(2016, 1, 1), datetime.datetime(2017, 12, 1)])
    plt.ylim([ymax - 2.5, ymax])
    plt.gcf().autofmt_xdate()
    plt.xlabel("DATE-OBS")
    plt.ylabel("Photometric Zeropoint %s" % select_filter)
    plt.title("Long term throughput  %s:%s in %s" % (select_site, select_telescope, select_filter))

    if instrument in ("kb97", "kb98"):
        plt.axvline(x=datetime.datetime(2017, 6, 30), color='k', linestyle='--')
    if instrument in ("fl03", "xxx"):
        plt.axvline(x=datetime.datetime(2017, 8, 31), color='k', linestyle='--')

    outfigname = "%s/photzptrend-%s-%s-%s.png" % (
        context.imagedbPrefix, select_site, select_telescope, select_filter)
    plt.savefig(outfigname, dpi=600)
    plt.close()

    plt.figure()
    plt.hist(zpsigselect, 50, range=[0, 1], normed=True)
    outerrorhistfname = "%s/errorhist-%s-%s-%s.png" % (
        context.imagedbPrefix, select_site, select_telescope, select_filter)
    plt.savefig(outerrorhistfname)
    plt.close()

    plt.figure()
    plt.plot(airmasselect, zpselect, ".", c="grey")
    plt.plot(airmasselect, zp_air, ".", c="blue")
    plt.xlabel("Airmass")
    plt.ylabel("Photomertic Zeropoint %s" % select_filter)
    meanzp = np.nanmedian(zpselect)
    plt.ylim([meanzp - 0.5, meanzp + 0.5])
    plt.savefig(
        "%s/airmasstrend-%s-%s-%s.png" % (context.imagedbPrefix, select_site, select_telescope, select_filter))
    plt.close()

    # Color terms
    plt.figure()
    selection = selection & np.logical_not(np.isnan(data['colorterm']))
    selection = selection & (np.abs(data['colorterm']) < 0.3)
    selection_lonoise = selection & (data['zpsig'] < 0.2)
    selection_hinoise = selection & (data['zpsig'] >= 0.2)

    plt.plot(data['dateobs'][selection_hinoise], data['colorterm'][
        selection_hinoise], '.', markersize=2, c="grey",
             label="color term [ hi sigma] %s " % select_filter)
    colortermselect = data['colorterm'][selection_lonoise]
    dateselect = data['dateobs'][selection_lonoise]
    meancolorterm = np.median(colortermselect)
    plt.plot(dateselect, colortermselect, 'o', markersize=2, c="blue",
             label="color term [low sigma] %s " % select_filter)
    plt.axhline(y=meancolorterm, color='r', linestyle='-')
    print("Color term filter %s : % 5.3f" % (select_filter, meancolorterm))

    if select_filter not in colorterms:
        colorterms[select_filter] = {}
    colorterms[select_filter][instrument] = meancolorterm

    plt.xlim([datetime.datetime(2016, 1, 1), datetime.datetime(2017, 12, 1)])
    plt.ylim([-0.2, 0.2])

    plt.savefig(
        "%s/colortermtrend-%s-%s-%s.png" % (context.imagedbPrefix, select_site, select_telescope, select_filter))
    plt.close()


def findUpperEnvelope(dateobs, datum, ymax=24.2):
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
        # Calculate the best throughput of a day
        todayzps = y[
            (x > startdate) & (x < startdate + datetime.timedelta(days=1)) & (
                y < ymax) & (y is not np.nan)]

        if len(todayzps) > 3:  # require a minimum amount of data for a night

            todayzps = np.sort(todayzps)[3:]
            maxzp = np.nanmax(todayzps)
            upperEnv = np.nanmean(todayzps[todayzps > (maxzp - stderror)])

            if upperEnv is not np.nan:
                day_x.append(startdate)
                day_y.append(upperEnv)

        startdate = startdate + datetime.timedelta(days=1)

    # filter the daily zero point variation. Work in progress.
    medianrange = 7
    newday_y = scipy.signal.medfilt(day_y, medianrange)

    return np.asarray(day_x), newday_y


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


def plotallmirrormodels(context, type='[2m0|1m0]'):
    import glob

    myfilter = 'rp'
    modellist = glob.glob("%s/mirrormodel*%s[abc]-%s.dat" % (context.imagedbPrefix, type, myfilter))

    for model in modellist:
        print(model)
        try:
            data = ascii.read(model, names=("date", "time", "zp"))
            datestring = np.core.defchararray.add(data['date'], "T")
            datestring = np.core.defchararray.add(datestring, data['time'])

            date = astt.Time(datestring, scale='utc', format='isot').to_datetime()
        except:
            continue

        plt.gcf().autofmt_xdate()
        plt.plot(date, data['zp'], label=model[-20:-8].replace('-', ':'))

    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', ncol=1)
    plt.xlabel('DATE-OBS')
    plt.ylabel("phot zeropoint %s" % myfilter)
    plt.xlim([datetime.datetime(2016, 1, 1), datetime.datetime(2017, 12, 1)])
    plt.title("Photometric zeropoint model in filter %s" % myfilter)
    plt.grid(True, which='both')
    plt.savefig("%s/allmodels_%s.png" % (context.imagedbPrefix, type), bbox_inches='tight')
    plt.close()


def parseCommandLine():
    """ Read command line parameters
    """

    parser = argparse.ArgumentParser(
        description='Calculate long-term trends in photometric database.')

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--databasedirectory', dest='imagedbPrefix', default='~/lcozpplots',
                        help='Directory containing photometryc databases')
    parser.add_argument('--site', dest='site', default=None, help='sites code for camera')
    parser.add_argument('--telescope', default=None,
                        help='Telescope id. written inform enclosure:telescope, e.g., "domb:1m0a"')
    parser.add_argument('--filter', default='rp', help='Which filter to process.', choices=['gp', 'rp', 'ip', 'zp'])

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    args.imagedbPrefix = os.path.expanduser(args.imagedbPrefix)

    return args


if __name__ == '__main__':
    plt.style.use('ggplot')
    matplotlib.rcParams['savefig.dpi'] = 600

    args = parseCommandLine()

    if args.site is not None:
        crawlsites = [args.site, ]
    else:
        crawlsites = telescopedict

    for site in crawlsites:

        if args.telescope is None:
            crawlScopes = telescopedict[site]
        else:
            crawlScopes = [args.telescope, ]

        for telescope in crawlScopes:
            print(site, telescope)
            plotlongtermtrend(site, telescope, args.filter, args, )

    plotallmirrormodels(args)
    plotallmirrormodels(args, type='0m4')

    sys.exit(0)
