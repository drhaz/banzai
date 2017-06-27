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

__author__ = 'dharbeck'

class PhotCalib(Stage):

    # These filters can be readilly transformed to the SDSS griz system with PS as a reference base.
    VALID_FILTERS = ('gp', 'rp', 'ip', 'zp')

    # To be replaced with map:
    # LCO filter -> sdss filter name, sdsss g-i color term, airmass term, default zero point.
    # possibly need this for all site:telescope:camera:filter settings. Maybe start with one
    # default and see where it goes.



    def __init__(self, pipeline_context):
        super(PhotCalib, self).__init__(pipeline_context)

    def inInCatalog (dec):
        return dec >= -30.0

    def do_stage(self, images):

        for i, image in enumerate(images):
            logging_tags = logs.image_config_to_tags(image, self.group_by_keywords)



class PS1Catalog:

    localChacheDir = None


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
                    maxsources=50000,
                    server=('http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx')):
        """
        Query Pan-STARRS DR1 @ MAST
        parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
                                          radius in degrees
                    mindet: minimum number of detection (optional)
                    maxsources: maximum number of sources
                 server: servername
        returns: astropy.table object

        copied from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
        # """
        r = requests.get(server,
        params= {'CAT':"PS1V3OBJECTS", 'RA': ra_deg, 'DEC': dec_deg,
             'SR': rad_deg, 'MAXOBJ': maxsources,
              'FORMAT': 'CSV',
              'ndetections': ('>%d' % mindet)})


        #write query data into local file
        #outf = open('panstarrs.dat', 'w')
        #outf.write(r.text)
        #outf.close()

        #input =  'panstarrs.dat'

        # need at least a few lines of input to reasoanably proceed.
        if (r.text.count('\n') > 10):
            input = StringIO.StringIO(r.text)
            table = np.genfromtxt (input,  names=True, skip_header=1, delimiter=',', )


            self.PS1toSDSS(table)

            return table
        return None



    def loadFitsCatalog (self, image):

        testimage = fits.open (image)

        ra  = testimage['SCI'].header['CRVAL1']
        dec = testimage['SCI'].header['CRVAL2']
        exptime = testimage['SCI'].header['EXPTIME']
        filterName = testimage['SCI'].header['FILTER']

        referenceInformation = self.FILTERMAPPING[filterName]

        if (referenceInformation is None):
            referenceInformation = self.FILTERMAPPING['rp']

        referenceFilterName = referenceInformation['refMag']


        instCatalog = testimage['CAT'].data
        image_wcs = WCS(testimage['SCI'].header)
        ras, decs = image_wcs.all_pix2world(instCatalog['x'], instCatalog['y'], 1)

        reftable = self.panstarrs_query(ra,dec,0.3)
        if reftable is None:
            print ("Failure on image %s, no reference table received." % (image))
            return
        # TODO: Check if catalog query was successful.

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



        return referenceRA, referenceDEC, instrumentalMag, referenceMag, referenceColor, matchDistance




# Example query

from astropy import units as u

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import re

def analyzeImage (imageName):


    outbasename = re.sub ('.fits.fz','', imageName)


    ps1 = PS1Catalog ()
    ra, dec, instmag, refmag, refcol, matchDist = ps1.loadFitsCatalog(imageName)



    magZP = refmag - instmag

    f=plt.figure()
    plt.plot (refmag, magZP, '.')
    plt.xlim([10,22])
    plt.ylim ([20,26])

    photzp = np.median (magZP)
    print ("Photometric zeropoint: %s  %5.2f" % (outbasename, photzp))

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

import glob

inputlist = glob.glob("../testing/lcogtdata-20170627-28/*.fits.fz")

for image in inputlist:


    analyzeImage(image)