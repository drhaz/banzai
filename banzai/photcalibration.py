from __future__ import absolute_import, division, print_function, unicode_literals
from banzai.stages import Stage
from banzai import logs
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
                    maxsources=5,
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
        r = requests.get(server,
        params= {'CAT':"PS1V3OBJECTS", 'RA': ra_deg, 'DEC': dec_deg,
             'SR': rad_deg, 'MAXOBJ': maxsources,
             'FORMAT': 'CSV',
             'ndetections': ('>%d' % mindet)})

        # write query data into local file

        print (r.text + "\n")
        table = np.genfromtxt (StringIO.StringIO(r.text),  names=True, skip_header=1, delimiter=',', )

        print (table['iMeanPSFMag'][0])
        self.PS1toSDSS(table)
        print (table['iMeanPSFMag'][0])
        return table

# Example query
ps1 = PS1Catalog ()
table = ps1.panstarrs_query(12.345, 67.89, 0.1)