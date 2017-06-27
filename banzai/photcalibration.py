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



    def __init__(self, pipeline_context):
        super(PhotCalib, self).__init__(pipeline_context)

    def inInCatalog (dec):
        return dec >= -30.0

    def do_stage(self, images):

        for i, image in enumerate(images):
            logging_tags = logs.image_config_to_tags(image, self.group_by_keywords)



class PS1Catalog:

    localChacheDir = None
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

        print ("converting PS1 to sdss system via finkbeiner")

        ps1colorterms = {}
        ps1colorterms['gMeanPSFMag'] = [-0.01808,-0.13595, 0.01941,-0.00183][::-1]
        ps1colorterms['rMeanPSFMag'] = [-0.01836,-0.03577, 0.02612,-0.00558][::-1]
        ps1colorterms['iMeanPSFMag'] = [ 0.01170,-0.00400, 0.00066,-0.00058][::-1]
        ps1colorterms['zMeanPSFMag'] = [-0.01062, 0.07529,-0.03592, 0.00890][::-1]
        pscolor = table['gMeanPSFMag'] - table['iMeanPSFMag']

        for filter in ps1colorterms:
            colorcorrection = np.polyval (ps1colorterms[filter], pscolor)
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
        """
        # r = requests.get(server,
        # params= {'CAT':"PS1V3OBJECTS", 'RA': ra_deg, 'DEC': dec_deg,
        #      'SR': rad_deg, 'MAXOBJ': maxsources,
        #      'FORMAT': 'CSV',
        #      'ndetections': ('>%d' % mindet)})
        #
        # input = StringIO.StringIO(r.text)
        # write query data into local file
        #outf = open('panstarrs.dat', 'w')
        #outf.write(r.text)
        #outf.close()

        input =  'panstarrs.dat'
        table = np.genfromtxt (input,  names=True, skip_header=1, delimiter=',', )


        self.PS1toSDSS(table)

        return table

# Example query

from astropy import units as u

def loadFitsCatalog (image):

    testimage = fits.open ('../testing/coj1m011-fl12-20170613-0073-e11.fits.fz')

    ra  = testimage['SCI'].header['CRVAL1']
    dec = testimage['SCI'].header['CRVAL2']
    exptime = testimage['SCI'].header['EXPTIME']
    phot = testimage['CAT'].data
    image_wcs = WCS(testimage['SCI'].header)
    ras, decs = image_wcs.all_pix2world(phot['x'], phot['y'], 1)


    return ra,dec,phot, ras, decs, exptime



ra,dec,phot,ras,decs, exptime = loadFitsCatalog   (None)
ps1 = PS1Catalog ()
reftable = ps1.panstarrs_query(ra,dec,0.3)



cphot = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
cref  = SkyCoord(ra=reftable['raMean']*u.degree, dec=reftable['decMean']*u.degree)
idx, d2d,d3d = cref.match_to_catalog_sky(cphot)
print (phot)
phott = np.transpose(phot)[idx]
phot  = np.transpose (phott)



distance = cref.separation (cphot[idx]).arcsecond

condition = (distance < 5) & (phot['FLUX'] > 0) & (reftable['rMeanPSFMag'] > 0)& (reftable['rMeanPSFMag'] < 26)

print ( distance )
instmag =  -2.5 * np.log10(phot['FLUX'] / exptime)


delta_r = ( reftable['rMeanPSFMag']  - instmag )[condition]
catalog_gi = (reftable['gMeanPSFMag']- reftable['iMeanPSFMag'])[condition]
catalog_r = reftable['rMeanPSFMag'][condition]

matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.plot (catalog_r, delta_r, '.')
plt.xlim([10,22])
plt.ylim ([20,26])

photzp = np.median (delta_r)
print ("Photometric zeropoint: %5.2f" % (photzp))

plt.axhline(y=photzp, color='r', linestyle='-')

plt.xlabel ("Reference catalog mag")
plt.ylabel ("Reference Mag - Instrumnetal Mag")
plt.title ("Photometric zeropoint %5.2f" % (photzp))
plt.savefig ("photoplot_zp.png")


plt.plot (catalog_gi, delta_r, '.')
plt.xlim([-0.5,3])
plt.ylim ([20,30])
plt.savefig ("photoplot_color.png")


