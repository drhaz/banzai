from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib

from banzai.stages import Stage
from banzai import logs
from banzai import images
from banzai.utils import image_utils
import os, subprocess, shlex
from astropy.io import fits
import tempfile
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units
import numpy as np
import requests
import StringIO
from astropy.io.votable import parse_single_table


import pickle

__author__ = 'dharbeck'

class PhotCalib():

    # These filters can be readilly transformed to the SDSS griz system with PS as a reference base.
    VALID_FILTERS = ('gp', 'rp', 'ip', 'zp')

    # To be replaced with map:
    # LCO filter -> sdss filter name, sdsss g-i color term, airmass term, default zero point.
    # possibly need this for all site:telescope:camera:filter settings. Maybe start with one
    # default and see where it goes.



  #  def __init__(self, pipeline_context):
        #super(PhotCalib, self).__init__(pipeline_context)

    def inInCatalog (dec):
        return dec >= -30.0

    def do_stage(self, images):

        for i, image in enumerate(images):
            logging_tags = logs.image_config_to_tags(image, self.group_by_keywords)



    def loadFitsCatalog (self, image):

        """ Load the photometry catalog from image 'CAT' extension and queries a PS1 matching catalog.
        If query is successful, return a set of matched catalog items


        :param image: input fits image path+name
        :return:
        """

        testimage = fits.open (image)

        # Boilerplate status infomration
        ra  = testimage['SCI'].header['CRVAL1']
        dec = testimage['SCI'].header['CRVAL2']
        exptime = testimage['SCI'].header['EXPTIME']
        filterName = testimage['SCI'].header['FILTER']
        airmass = testimage['SCI'].header['AIRMASS']

        if filterName not in PS1Catalog.FILTERMAPPING:
            print ("Filter not viable for photometrric calibration. Sorry")
            return

        referenceInformation = PS1Catalog.FILTERMAPPING[filterName]

        if (referenceInformation is None):
            referenceInformation = PS1Catalog.FILTERMAPPING['rp']

        referenceFilterName = referenceInformation['refMag']


        # Load photometry catalog from image, and transform into RA/Dec coordinates
        instCatalog = testimage['CAT'].data
        image_wcs = WCS(testimage['SCI'].header)
        ras, decs = image_wcs.all_pix2world(instCatalog['x'], instCatalog['y'], 1)

        # Query reference catalog
        pscatalog = PS1Catalog()
        reftable = pscatalog.panstarrs_query(ra,dec,0.3)

        if reftable is None:
            print ("Failure on image %s, no reference table received." % (image))
            return

        # Start the catalog matching, using astropy skycoords built-in functions.
        cInstrument = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
        cReference  = SkyCoord(ra=reftable['raMean']*u.degree, dec=reftable['decMean']*u.degree)
        idx, d2d,d3d = cReference.match_to_catalog_sky(cInstrument)

        # reshuffle the source catalog to index-match the reference catalog.
        #  There is probably as smarter way of doing this!
        instCatalogT = np.transpose(instCatalog)[idx]
        instCatalog  = np.transpose (instCatalogT)

        # Measure the distance between matched pairs. Important to down-select viable pairs.
        distance = cReference.separation (cInstrument[idx]).arcsecond


        # Define a reasonable condition on what is a good match on good photometry
        condition = (distance < 5) & (instCatalog['FLUX'] > 0) & (reftable[referenceFilterName] > 0)& (reftable[referenceFilterName] < 26)
        if (exptime <=0):
            extime =1
        # Calculate instrumental magnitude from PSF instrument photometry
        instmag =  -2.5 * np.log10(instCatalog['FLUX'][condition] / exptime)

        # Calculate the magnitude difference between reference and inst catalog
        instrumentalMag =  instmag
        referenceColor  = (reftable['gMeanPSFMag']- reftable['iMeanPSFMag'])[condition]

        referenceMag = reftable[referenceFilterName][condition]
        referenceRA  = reftable['raMean'][condition]
        referenceDEC = reftable['decMean'][condition]
        matchDistance = distance[condition]


        return referenceRA, referenceDEC, instrumentalMag, referenceMag, referenceColor, matchDistance, referenceFilterName, filterName, airmass


    def analyzeImage (self, imageName, pickle="photzp.db"):
        """ Do full photometric zeropoint analysis on an image"""

        outbasename = re.sub ('.fits.fz','', imageName)



        ra, dec, instmag, refmag, refcol, matchDist, refFilter, instFilter, airmass = self.loadFitsCatalog(imageName)

        magZP = refmag - instmag

        f=plt.figure()
        plt.plot (refmag, magZP, '.')
        plt.xlim([10,22])
        plt.ylim ([20,26])

        photzp = np.median (magZP)
        print ("Photometric zeropoint: %s %s %5.2f %5.2f" % (outbasename, instFilter, airmass, photzp))



        plt.axhline(y=photzp, color='r', linestyle='-')

        plt.xlabel ("Reference catalog mag")
        plt.ylabel ("Reference Mag - Instrumnetal Mag")
        plt.title ("Photometric zeropoint %s %5.2f" % (outbasename, photzp))
        plt.savefig ("%s_zp.png" % (outbasename))
        plt.close()

        f=plt.figure()
        plt.plot (refcol , magZP - photzp, '.')
        plt.xlim([-0.5,3])
        plt.ylim ([-1,1])
        plt.xlabel ("(g-r)_{SDSS} Reference")
        plt.ylabel ("Reference Mag - Instrumnetal Mag - ZP (%5.2f)" %( photzp))
        plt.title ("Color correction %s " % (outbasename))
        plt.savefig ("%s_color.png" % (outbasename))
        plt.close()

        f=plt.figure()
        plt.scatter (ra,dec, c= magZP - photzp, vmin=-0.2, vmax=0.2, edgecolor='none',
                     s=9, cmap=matplotlib.pyplot.cm.get_cmap( 'nipy_spectral') )
        plt.colorbar()
        plt.title ("Spacial variation of phot. Zeropoint %s" % (outbasename))
        plt.xlabel ("RA")
        plt.ylabel ("Dec")
        plt.savefig ("%s_zpmap.png" % (outbasename))
        plt.close()


class PS1Catalog:

    localCacheDir = None

    FILTERMAPPING = {}
    FILTERMAPPING['gp'] = {'refMag': 'gMeanPSFMag', 'colorTerm': 0.0, 'airmassTerm': 0.0, 'defaultZP': 0.0}
    FILTERMAPPING['rp'] = {'refMag': 'rMeanPSFMag', 'colorTerm': 0.0, 'airmassTerm': 0.0, 'defaultZP': 0.0}
    FILTERMAPPING['ip'] = {'refMag': 'iMeanPSFMag', 'colorTerm': 0.0, 'airmassTerm': 0.0, 'defaultZP': 0.0}
    FILTERMAPPING['zp'] = {'refMag': 'zMeanPSFMag', 'colorTerm': 0.0, 'airmassTerm': 0.0, 'defaultZP': 0.0}

    ps1colorterms = {}
    ps1colorterms['gMeanPSFMag'] = [-0.01808,-0.13595, 0.01941,-0.00183][::-1]
    ps1colorterms['rMeanPSFMag'] = [-0.01836,-0.03577, 0.02612,-0.00558][::-1]
    ps1colorterms['iMeanPSFMag'] = [ 0.01170,-0.00400, 0.00066,-0.00058][::-1]
    ps1colorterms['zMeanPSFMag'] = [-0.01062, 0.07529,-0.03592, 0.00890][::-1]


    def __init__ (self):
        pass

    def PS1toSDSS (self, table):
        """ PS1 catalog is calibrated to PS12 photometric system, which is different from SDSS

        This procedure transforms a catalog from the PS1 system into the SDSS catalog following
        Finkbeiner 2016
        http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/meta#apj522061s2-4 Table 2

        Note that this transformation is valid for stars only. For the purpose of photometric
        calibration, it is desirable to select point sources onyl from the input catalog.
        """
        pscolor = table['gMeanPSFMag'] - table['iMeanPSFMag']

        for filter in self.ps1colorterms:
            colorcorrection = np.polyval (self.ps1colorterms[filter], pscolor)
            table[filter] -= colorcorrection
        return table


    def panstarrs_query(self, ra_deg, dec_deg, rad_deg, mindet=5,
                    maxsources=5000,
                    server=('http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx')):

        """ Online query to PS1 from MAST

            Returns a table as it is retrieved from MAST, loaded into numpy.genfromtxt().
            The PS photometry values xMeanSPFMag are converted into the SDSSS photometric system.

        """
        print ("Start retrieve catalog data at ra % 8.4f dec % 8.4f" % (ra_deg,dec_deg))
        r = requests.get(server,
        params= {'CAT':"PS1V3OBJECTS", 'RA': ra_deg, 'DEC': dec_deg,
             'SR': rad_deg, 'MAXOBJ': maxsources,
              'FORMAT': 'CSV',
              'ndetections': ('>%d' % mindet)})

        # need at least a few lines of input to reasonably proceed.

        nlines = r.text.count('\n')
        if (nlines > 10):
            input = StringIO.StringIO(r.text)
            table = np.genfromtxt (input,  names=True, skip_header=1, delimiter=',', )

            self.PS1toSDSS(table)

            print ("\tDone retrieve table, %d lines" % nlines)
            return table
        print ("\t failed retrieve table)")
        return None


# Example query

from astropy import units as u

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import re



import glob

inputlist = glob.glob("../testing/lcogtdata-20170627-28/*.fits.fz")


photzpStage = PhotCalib()


for image in inputlist:


    photzpStage.analyzeImage(image)