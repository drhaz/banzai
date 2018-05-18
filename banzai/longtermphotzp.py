from __future__ import absolute_import, division, print_function, unicode_literals
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import datetime
import sys
import calendar
import astropy.time as astt
import scipy.signal
import argparse
import logging
import glob
import os
import sqlite3
from astropy.io import ascii
from astropy.table import Table
from itertools import cycle

assert sys.version_info >= (3,5)
_logger = logging.getLogger(__name__)



airmasscorrection = {'gp': 0.17, 'rp': 0.09, 'ip': 0.06, 'zp': 0.05, }

starttime = datetime.datetime(2016, 1, 1)
#starttime = datetime.datetime(2018,3,1)

endtime   = datetime.datetime.utcnow().replace(day=28) + datetime.timedelta(days=31+4)
endtime.replace(day =1)



colorterms = {}
telescopedict = {
    'lsc': ['doma-1m0a', 'domb-1m0a', 'domc-1m0a', 'aqwa-0m4a', 'aqwb-0m4a'],
    'coj': ['clma-2m0a', 'doma-1m0a', 'domb-1m0a', 'clma-0m4a', 'clma-0m4b'],
    'ogg': ['clma-2m0a', 'clma-0m4b', 'clma-0m4c'],
    'elp': ['doma-1m0a', 'aqwa-0m4a'],
    'cpt': ['doma-1m0a', 'domb-1m0a', 'domc-1m0a', 'aqwa-0m4a'],
    'tfn': ['aqwa-0m4a', 'aqwa-0m4b'],
    'sqa': ['doma-0m8a'],
    'bpl': ['doma-1m0a']
}


telescopecleaning = {
    'lsc-doma-1m0a' : [datetime.datetime(2018, 4, 5),] ,
    'lsc-domb-1m0a' : [datetime.datetime(2018, 4, 5),] ,
    'lsc-domc-1m0a' : [datetime.datetime(2017, 8, 31), datetime.datetime(2018, 4, 5),] ,
    'lsc-aqwa-0m4a' : [datetime.datetime(2018, 4, 17),] ,
    'lsc-aqwb-0m4a' : [datetime.datetime(2018, 4, 17),] ,
    'coj-clma-0m4a' : [datetime.datetime(2017, 6, 30),] ,
    'coj-clma-0m4b' : [datetime.datetime(2017, 6, 30),] ,
    'elp-doma-1m0a' : [datetime.datetime(2017, 9, 20), datetime.datetime(2018, 4, 5),datetime.datetime(2018, 5, 6), ] ,
    'ogg-clma-2m0a' : [datetime.datetime(2017, 10,20),],
    'cpt-doma-1m0a' : [datetime.datetime(2017, 11, 15),] ,
    'cpt-domb-1m0a' : [datetime.datetime(2017, 11, 15),] ,
    'cpt-domc-1m0a' : [datetime.datetime(2017, 11, 15),] ,
}


class photdbinterface:

    createstatement = "CREATE TABLE IF NOT EXISTS lcophot (" \
                      "name TEXT PRIMARY KEY, " \
                      "dateobs text," \
                      " site text," \
                      " dome text," \
                      " telescope text," \
                      " camera text," \
                      " filter text," \
                      " airmass real," \
                      " zp real," \
                      " colorterm real," \
                      " zpsig real)"



    def __init__(self, fname):
        _logger.info ("Open data base file %s" % (fname))
        self.sqlite_file = fname
        self.conn = sqlite3.connect(self.sqlite_file)
        self.conn.execute(self.createstatement)
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.commit()

    def addphotzp (self, datablob, commit = True) :

        _logger.info ("About to insert: %s" % str(datablob))

        with self.conn:
            self.conn.execute ("insert or replace into lcophot values (?,?,?,?,?,?,?,?,?,?,?)", datablob)

            if (commit):
                self.conn.commit()


    def exists(self, fname):
        """ Check if entry as identified by file name already exists in database
        """

        cursor = self.conn.execute ("select * from lcophot where name=? limit 1", (fname,))
        res = cursor.fetchone()
        if res is None:
            return False
        return (len (res) > 0)


    def readoldfile (self, oldname):
        """ Ingest a legacy ascii photomerty zeropoint file."""
        _logger.info ("Ingesting old style ascii datgabase : %s" % (oldname))
        data = readDataFile(oldname, False)

        for line in data:
            print (line)
            self.addphotzp(line, commit = False)
        self.conn.commit()


    def close(self):
        """ Clode the database safely"""

        _logger.info ("Closing data base file %s " % (self.sqlite_file))
        self.conn.commit()
        self.conn.close()


    def readRecords (self, site = None, dome = None, telescope = None, camera = None):
        """  Read the photometry records fromt eh database, optionally filtered by site, dome, telescope, and camera.

        """


        query = "select name,dateobs,site,dome,telescope,camera,filter,airmass,zp,colorterm,zpsig from lcophot " \
                "where (site like ? AND dome like ? AND telescope like ? AND camera like ?)"

        args = (site if site is not None else '%',
                dome if dome is not None else '%',
                telescope if telescope is not None else '%',
                camera if camera is not None else '%',)

        cursor = self.conn.execute(query, args)

        allrows = np.asarray(cursor.fetchall())
        if len (allrows) == 0:
            return None

        t = Table (allrows, names = ['name','dateobs','site','dome','telescope','camera','filter','airmass','zp','colorterm','zpsig'])
        t['dateobs'] = t['dateobs'].astype (str)
        t['dateobs'] = astt.Time(t['dateobs'], scale='utc', format=None).to_datetime()

        t['zp'] = t['zp'].astype(float)

        t['airmass'] [ 'UNKNOWN' == t['airmass'] ]  = 'nan'
        t['airmass'] = t['airmass'].astype(float)

        t['zpsig'] = t['zpsig'].astype(float)
        t['colorterm'] = t['colorterm'].astype(float)


        if 'fl06' in t['camera']:
            # fl06 was misconfigured with a wrong gain, which trickles down through the banzai processing.
            # The correct gain was validated Nov 27th 2017 on existing data.
            dateselect = ( t['dateobs'] < datetime.datetime(year=2017,month=11,day=17) ) & (t['camera'] == 'fl06')
            t['zp'][dateselect] = t['zp'][dateselect] - 2.5 * math.log10 (1.82 / 2.45)

        if 'fl05' in t['camera']:
        # fl06 was misconfigured with a wrong gain, which trickles down through the banzai processing.
        # The correct gain was validated Nov 27th 2017 on existing data.
            dateselect =  (t['camera'] == 'fl05')
            t['zp'][dateselect] = t['zp'][dateselect] - 2.5 * math.log10 (1.69 / 2.09)

        if 'fl11' in  t['camera']:
            #
            dateselect =  (t['camera'] == 'fl11')
            t['zp'][dateselect] = t['zp'][dateselect] - 2.5 * math.log10 (1.85 / 2.16)

        if 'kb96' in t['camera']:
            dateselect = ( t['dateobs'] > datetime.datetime(year=2017,month=11,day=15) ) & ( t['dateobs'] < datetime.datetime(year=2018,month=4,day=10) ) & (t['camera'] == 'kb96')
            t['zp'][dateselect] = t['zp'][dateselect] - 2.5 * math.log10 (0.851 / 2.74)

        return t


def readDataFile(inputfile, gaincorrect = True):
    print ("DEPRECATED!")
    exit (0)
    with open(inputfile, 'r') as file:
        contents = file.read()
        file.close()

        contents = contents.replace('UNKNOWN', 'nan')
        data = ascii.read(contents, names=['name', 'dateobs', 'site', 'dome',
                                           'telescope', 'camera', 'filter', 'airmass',
                                           'zp', 'colorterm', 'zpsig'], )

        data['dateobs'] = astt.Time(data['dateobs'], scale='utc', format='isot').to_datetime()

        if not gaincorrect:
            return data

        if 'fl06' in inputfile:
            # fl06 was misconfigured with a wrong gain, which trickles down through the banzai processing.
            # The correct gain was validated Nov 27th 2017 on existing data.
            dateselect = data['dateobs'] < datetime.datetime(year=2017,month=11,day=17)
            data['zp'][dateselect] = data['zp'][dateselect] - 2.5 * math.log10 (1.82 / 2.45)

        if 'fl11' in inputfile:
            # fl06 was misconfigured with a wrong gain, which trickles down through the banzai processing.
            # The correct gain was validated Nov 27th 2017 on existing data.
            data['zp'] = data['zp'] - 2.5 * math.log10 (1.85 / 2.16)

        return data

    return None


def getCombineddataByTelescope(site, telescope, context, instrument=None, cacheddb=None):
    """
    Concatenate all zeropoint data for a site, and select by telescope and instrument.
    :param site:
    :param telescope: string slecting dome *& telescope: 'domb:1m0a'
    :param context:
    :param instrument:
    :return: concatenated data for a site / tel / isntrument selection.
    """

    if cacheddb is None:
        db = photdbinterface(context.database)
    else:
        db=cacheddb

    print (site,telescope,instrument)
    dome, tel = telescope.split ("-")

    results =  db.readRecords(site,dome,tel,instrument)
    if cacheddb is None:
        db.close()
    return results


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
                print (site, telescope, instrument, data)

    if alldata is None:
        return None

    # Helpful diagnostic tool to see what data are pulled in.
    # domedict = np.unique (alldata['dome'])
    # for dome in domedict:
    #     teldict = np.unique(alldata [ alldata['dome'] == dome] ['telescope'])
    #     print (site, dome, teldict)

    dome, tel = telescope.split('-')
    selection = (alldata['dome'] == dome) & (alldata['telescope'] == tel)
    alldata = alldata[selection]
    return alldata


def plotlongtermtrend(select_site, select_telescope, select_filter, context, instrument=None, cacheddb = None):

    data = getCombineddataByTelescope(select_site, select_telescope, context, instrument, cacheddb=cacheddb)

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

    for telid in telescopecleaning:
        _site,_enc,_tel = telid.split ("-")

        if (_site == select_site) and (select_telescope == '%s-%s' % (_enc,_tel)):
            for event in telescopecleaning[telid]:
                plt.axvline (x=event, color='grey', linestyle='--')



    uniquecameras = np.unique(cameraselect)



    for uc in uniquecameras:
        # plot zeropoint with differnt markers per camera

        plt.plot(dateselect[(zpsigselect <= photzpmaxnoise) & (cameraselect == uc)],
                 zp_air[(zpsigselect <= photzpmaxnoise) & (cameraselect == uc)],
                 'o', markersize=2, label=uc)
        plt.plot(dateselect[zpsigselect > photzpmaxnoise], zp_air[zpsigselect > photzpmaxnoise], '.',
                 markersize=1, c="grey", label='rejected' )

    if _x is not None:
        plt.plot(_x, _y, "-", c='red', label='upper envelope')

    else:
        print("Mirror model failed to compute. not plotting !")




    plt.legend()
    plt.xlim([starttime, endtime])
    plt.ylim([ymax - 3.5, ymax])
    plt.gcf().autofmt_xdate()
    plt.xlabel("DATE-OBS")
    plt.ylabel("Photometric Zeropoint %s" % select_filter)
    plt.title("Long term throughput  %s:%s in %s" % (select_site, select_telescope, select_filter))


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
    plt.title ("Global airmass trend and correction check")

    meanzp = np.nanmedian(zpselect)
    plt.ylim([meanzp - 0.5, meanzp + 0.5])

    plt.savefig("%s/airmasstrend-%s-%s-%s.png" % (context.imagedbPrefix, select_site, select_telescope, select_filter))
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

    plt.xlim([starttime, endtime])
    plt.ylim([-0.2, 0.2])
    plt.gcf().autofmt_xdate()

    plt.title("Color term (g-r)  %s:%s in %s" % (select_site, select_telescope, select_filter))
    plt.xlabel("DATE-OBS")
    plt.ylabel("Color term coefficient (g-r)")
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

    stderror = 0.03

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

            todayzps = np.sort(todayzps)[1:]
            maxzp = np.nanmax(todayzps)
            upperEnv = np.nanmean(todayzps[todayzps > (maxzp - stderror)])

            if upperEnv is not np.nan:
                day_x.append(startdate)
                day_y.append(upperEnv)

        startdate = startdate + datetime.timedelta(days=1)

    # filter the daily zero point variation. Work in progress.
    medianrange = 9
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


def plotallmirrormodels(context, type='[2m0a|1m0a]', range=[22.5,25.5]):
    import glob

    myfilter = context.filter
    modellist = glob.glob("%s/mirrormodel*%s[abc]-%s.dat" % (context.imagedbPrefix, type, myfilter))
    print (modellist[0][16:-12])

    modellist.sort(key = lambda x: x[-20:-17] + x[-11:-7].replace("2","0")+x[-16:-12])
    print (modellist)

    plt.rc('lines', linewidth=1)
    prop_cycle=  cycle( ['-', '-.'])

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
        plt.plot(date, data['zp'],next(prop_cycle),  label=model[-20:-7].replace('-', ':'), )

    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', ncol=1)
    plt.xlabel('DATE-OBS')
    plt.ylabel("phot zeropoint %s" % myfilter)
    plt.xlim([starttime, endtime])
    plt.ylim(range)
    plt.title("Photometric zeropoint model in filter %s" % myfilter)
    plt.grid(True, which='both')
    plt.savefig("%s/allmodels_%s_%s.png" % (context.imagedbPrefix, type, context.filter), bbox_inches='tight')
    plt.close()


def parseCommandLine():
    """ Read command line parameters
    """

    parser = argparse.ArgumentParser(
        description='Calculate long-term trends in photometric database.')

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--outputdirectory', dest='imagedbPrefix', default='~/lcozpplots',
                        help='Directory containing photometryc databases')
    parser.add_argument('--database', default = '~/lcozpplots/lcophotzp.db')
    parser.add_argument('--site', dest='site', default=None, help='sites code for camera')
    parser.add_argument('--telescope', default=None,
                        help='Telescope id. written inform enclosure-telescope, e.g., "domb-1m0a"')
    parser.add_argument('--filter', default='rp', help='Which filter to process.', choices=['gp', 'rp', 'ip', 'zp'])

    parser.add_argument('--testdb',  action='store_true')
    parser.add_argument('--importold',  action='store_true')


    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    args.imagedbPrefix = os.path.expanduser(args.imagedbPrefix)
    args.database = os.path.expanduser(args.database)


    return args


import webbrowser
def renderHTMLPage (args):

    outputfile = "%s/index.html" % (args.imagedbPrefix)

    message = """<html>
<head></head>
<body><title>LCO Zeropoint Plots</title>

<h1> Overview </h1>
     <a href="allmodels_[2m0a|1m0a]_rp.png"> <img src="allmodels_[2m0a|1m0a]_rp.png"/ width="800"> </a>
     <a href="allmodels_0m4_rp.png"><img src="allmodels_0m4_rp.png" width="800"/></a>
    <p/>
    
<h1> Details by Site: </h1>
"""

    for site in telescopedict:
        message = message + " <h2> %s </h2>\n" % (site)

        zptrendimages = glob.glob ("%s/photzptrend-%s-????-????-rp.png" % (args.imagedbPrefix, site))

        zptrendimages.sort(key = lambda x: x[-16: -4])

        print (zptrendimages)
        for zptrend in zptrendimages:
            zptrend = zptrend.replace("%s/" % args.imagedbPrefix, "")
            line = '<a href="%s"><img src="%s" width="600"/></a>  <img src="%s" width="600"/>  <img src="%s" width="600"/><p/>' % (zptrend, zptrend, zptrend.replace('photzptrend', 'colortermtrend'), zptrend.replace('photzptrend', 'airmasstrend'))
            message = message + line


    message = message + "<p/>Figures updated %s UTC <p/></body></html>" % (datetime.datetime.utcnow())

    with open (outputfile, 'w+') as f:
        f.write (message)
        f.close()
        #webbrowser.open('file://%s' % outputfile, new=False)



if __name__ == '__main__':
    plt.style.use('ggplot')
    matplotlib.rcParams['savefig.dpi'] = 600

    args = parseCommandLine()


    if (args.testdb):
        db = photdbinterface(args.database)
        db.readRecords()
        db.close()
        exit()


    if (args.importold):
        print ("testing db interface")
        db = photdbinterface(args.database)
        dbfiles = glob.glob ("/home/dharbeck/lcozpplots/???-*.db")
        for dbfile in dbfiles:

            print ("Importing  " + dbfile)
            db.readoldfile(dbfile)

        db.close()
        exit(0)

    if args.site is not None:
        crawlsites = [args.site, ]
    else:
        crawlsites = telescopedict



    db = photdbinterface(args.database)


    if True:
       for site in crawlsites:
            if args.telescope is None:
                crawlScopes = telescopedict[site]
            else:
                crawlScopes = [args.telescope, ]

            for telescope in crawlScopes:
                print(site, telescope)

                plotlongtermtrend(site, telescope, args.filter, args, cacheddb=db )


    db.close()

    plotallmirrormodels(args)
    plotallmirrormodels(args, type='0m4', range=[20,23])


    # Make a fancy HTML page
    renderHTMLPage(args)

    sys.exit(0)

