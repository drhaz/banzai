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
from numpy import rec
from astropy.io.votable import parse_single_table

import logging
import pickle

__author__ = 'dharbeck'

class PhotCalib():

    # These filters can be readilly transformed to the SDSS griz system with PS as a reference base.
    VALID_FILTERS = ('gp', 'rp', 'ip', 'zp')

    # To be replaced with map:
    # LCO filter -> sdss filter name, sdsss g-i color term, airmass term, default zero point.
    # possibly need this for all site:telescope:camera:filter settings. Maybe start with one
    # default and see where it goes.
    ps1 = None


    def __init__(self):
        #super(PhotCalib, self).__init__(pipeline_context)
        self.ps1 = PS1Catalog("/home/dharbeck/Catalogs/PS1")

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
        dateobs = testimage['SCI'].header['DATE-OBS']

        if filterName not in self.ps1.FILTERMAPPING:
            print ("Filter not viable for photometrric calibration. Sorry")
            return

        referenceInformation = self.ps1.FILTERMAPPING[filterName]

        if (referenceInformation is None):
            referenceInformation = self.ps1.FILTERMAPPING['rp']

        referenceFilterName = referenceInformation['refMag']


        # Load photometry catalog from image, and transform into RA/Dec coordinates
        instCatalog = testimage['CAT'].data
        image_wcs = WCS(testimage['SCI'].header)
        ras, decs = image_wcs.all_pix2world(instCatalog['x'], instCatalog['y'], 1)

        # Query reference catalog

        reftable = self.ps1.panstarrs_query(ra-0.2,dec-0.2,ra+0.2,dec+0.2)

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


        return referenceRA, referenceDEC, instrumentalMag, referenceMag, referenceColor, matchDistance, referenceFilterName, filterName, airmass, dateobs


    def analyzeImage (self, imageName, pickle="photzp.db"):
        """ Do full photometric zeropoint analysis on an image"""

        outbasename = re.sub ('.fits.fz','', imageName)



        ra, dec, instmag, refmag, refcol, matchDist, refFilter, instFilter, airmass, dateobs = self.loadFitsCatalog(imageName)

        magZP = refmag - instmag

        f=plt.figure()
        plt.plot (refmag, magZP, '.')
        plt.xlim([10,22])
        plt.ylim ([20,26])

        photzp = np.median (magZP)




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

        with open(pickle,'a') as f:
            output = "%s %s %s %s % 6.3f \n" % (imageName, dateobs, airmass, instFilter, photzp)
            print (output)
            f.write(output)
            f.close()
        return photzp

import os

import math
class PS1IPP:

    def get_reference_catalog(self, ra, dec, radius, basedir, overwrite_select=False ):

        logger = logging.getLogger("ReadCatalog")
        catalog_filenames = None

        # print "In get_ref_catalog, cattype=%s, dir=%s" % (cattype, basedir)

        if (basedir is None or not os.path.isdir(basedir)):
            #(basedir is not None and not os.path.isdir(basedir))):
            logger.warning("Unable to find reference catalog: %s" % (str(basedir)))
            return None

        # Load the SkyTable so we know in what files to look for the catalog"
        logger.debug("Using catalog found in %s" % (basedir))
        skytable_filename = "%s/SkyTable.fits" % (basedir)
        if (not os.path.isfile(skytable_filename)):
            logger.error("Unable to find catalog index file in %s!" % (basedir))
            return None

        skytable_hdu = fits.open(skytable_filename)

        select_header_list = None
        if ('NSELECT' in skytable_hdu[0].header):
            n_select = skytable_hdu[0].header['NSELECT']
            select_header_list = [None] * n_select
            for i in range(n_select):
                select_header_list[i] = skytable_hdu[0].header['SELECT%02d' % (i+1)]
            logger.debug("Selecting the following columns: %s" % (", ".join(select_header_list)))

        #print skytable_hdu.info()

        skytable = skytable_hdu['SKY_REGION'].data
        #print skytable[:3]

        # Select entries that match our list
        # print ra, dec, radius, type(ra), type(dec), type(radius)
        #logger.debug("# Searching for stars within %.1f degress around %f , %f ..." % (radius, ra, dec))

        if (not radius == None and radius > 0):
            min_dec = dec - radius
            max_dec = dec + radius
            min_ra = ra - radius/math.cos(math.radians(dec))
            max_ra = ra + radius/math.cos(math.radians(dec))
        else:
            min_dec, max_dec = dec[0], dec[1]
            min_ra, max_ra = ra[0], ra[1]

        logger.debug("Querying catalog: Ra=%f...%f Dec=%f...%f" % (min_ra, max_ra, min_dec, max_dec))

        if (max_ra > 360.):
            # This wraps around the high end, shift all ra values by -180
            # Now all search RAs are ok and around the 180, next also move the catalog values
            selected = skytable['R_MIN'] < 180
            skytable['R_MAX'][selected] += 360
            skytable['R_MIN'][selected] += 360
        if (min_ra < 0):
            # Wrap around at the low end
            selected = skytable['R_MAX'] > 180
            skytable['R_MAX'][selected] -= 360
            skytable['R_MIN'][selected] -= 360

        if (True): print ("# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec))

        try:
            needed_catalogs = (skytable['PARENT'] > 0) & (skytable['PARENT'] < 25) & \
                           (skytable['R_MAX'] > min_ra) & (skytable['R_MIN'] < max_ra) & \
                            (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)
        except KeyError:
            # try without the PARENT field
            needed_catalogs =   (skytable['R_MAX'] > min_ra)  & (skytable['R_MIN'] < max_ra) & \
                            (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)

        #print skytable[needed_catalogs]

        files_to_read = skytable['NAME'][needed_catalogs]
        files_to_read = [f.strip() for f in files_to_read]
        logger.debug(files_to_read)

        # Now quickly go over the list and take care of all filenames that still have a 0x00 in them
        for i in range(len(files_to_read)):
            found_at = files_to_read[i].find('\0')
            if (found_at > 0):
               files_to_read[i] = files_to_read[i][:found_at]

        # Now we are with the skytable catalog, so close it
        skytable_hdu.close()
        del skytable

        #print files_to_read

        # Load all frames, one by one, and select all stars in the valid range.
        # Then add them to the catalog with RAs and DECs
        full_catalog = None #numpy.zeros(shape=(0,6))
        catalog_filenames = []

        # Start iterating though catalogs
        for catalogname in files_to_read:

            catalogfile = "%s/%s" % (basedir, catalogname)
            # print catalogfile
            if (not os.path.isfile(catalogfile)):
                # not a file, try adding .fits to the end of the filename
                if (os.path.isfile(catalogfile+".fits")):
                    catalogfile += ".fits"
                else:
                    # neither option (w/ or w/o .fits added is a file)
                    logger.warning("Catalog file (%s) not found (base-dir: %s)" % (os.path.abspath(catalogfile), basedir))
                    continue

            try:
                hdu_cat = fits.open(catalogfile)
            except:
                logger.warning("Unable to open catalog file %s" % (catalogfile))
                continue

            catalog_filenames.append(catalogfile)
            logger.debug("Adding %s to list of catalog files being used" % (catalogfile))

            # read table into a nd-array buffer
            cat_full = self.table_to_ndarray(hdu_cat[1], select_header_list=select_header_list,
                                        overwrite_select=overwrite_select)
            # print cat_full.shape

            # Read the RA and DEC values
            cat_ra  = cat_full[:,0]
            cat_dec = cat_full[:,1]

            # To slect the right region, shift a temporary catalog
            cat_ra_shifted = cat_ra
            if (max_ra > 360.):
                cat_ra_shifted[cat_ra < 180] += 360
            elif (min_ra < 0):
                cat_ra_shifted[cat_ra > 180] -= 360

            select_from_cat = (cat_ra_shifted > min_ra) & (cat_ra_shifted < max_ra ) & (cat_dec > min_dec) & (cat_dec < max_dec)

            array_to_add = cat_full[select_from_cat]
            logger.debug("Read %d sources from %s" % (array_to_add.shape[0], catalogname))


            if (full_catalog is None):
                full_catalog = array_to_add
            else:
                full_catalog = np.append(full_catalog, array_to_add, axis=0)
                #print photom_grizy[:3,:]

            if (full_catalog is None):
                logger.warning("No stars found in area %s, %s from catalog %s" % (
                    str(ra), str(dec),
                    #ra[0], ra[1], dec[0], dec[1],
                    basedir))
            else:
                logger.debug("Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read)))




        return full_catalog

    def table_to_ndarray(self, tbhdu, select_header_list=None, overwrite_select=False):

        logger = logging.getLogger("Table2Array")
        n_entries = tbhdu.header['NAXIS2']

        if (select_header_list is None):

            # no selection of what columns to read

            n_fields = tbhdu.header['TFIELDS']

            logger.debug("Found %d columns and %d rows" % (n_fields, n_entries))

            databuffer = np.empty((n_entries, n_fields))
            for i in range(n_fields):
                try:
                    coldata = np.array(tbhdu.data.field(i))
                    databuffer[:,i] = coldata[:] #tbhdu.data.field(i)
                except ValueError:
                    pass

        else:

            n_fields = len(select_header_list)
            if (overwrite_select):
                n_fields += tbhdu.header['TFIELDS']
                logger.debug("Reading select list of %d columns for %d sources" % (n_fields, n_entries))

            databuffer = np.empty((n_entries, n_fields))
            for i, fieldname in enumerate(select_header_list):
                try:
                    coldata = np.array(tbhdu.data.field(fieldname))
                    databuffer[:, i] = coldata[:]  # tbhdu.data.field(i)
                except ValueError:
                    pass

            if (overwrite_select):
                for i in range(tbhdu.header['TFIELDS']):
                    try:
                        coldata = np.array(tbhdu.data.field(i))
                        databuffer[:,i+len(select_header_list)] = coldata[:] #tbhdu.data.field(i)
                    except ValueError:
                        pass


        # print databuffer.shape
        # numpy.savetxt(fn[:-5]+".dump", databuffer)
        return databuffer


class PS1Catalog:

    localCacheIndexBaseDir = "/home/dharbeck/Catalogs/PS1"
    indexTable = None
    gridsize=0.5

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


    def __init__ (self, cacheIndex = None):

        if cacheIndex is not None:
            self.localCacheIndexBaseDir = cacheIndex
            print ("PS1 cache detected: %s" % self.localCacheIndexBaseDir)

            data = []

            for dec in np.arange (-30,+90,self.gridsize):
                idx = 0
                for ra in np.arange (0,360,self.gridsize):
                    name = "%04d/%05d-%010d.txt" % (dec %10, dec*100, idx)
                    idx = idx + 1
                    data.append ( (ra, ra+self.gridsize, dec, dec+self.gridsize, str(name)))

            dtype = np.dtype({'names':['R_MIN','R_MAX','D_MIN','D_MAX','NAME'], 'formats':[float,float,float,float,'S25']})

            self.indexTable = np.array(data,  dtype=dtype )

            pass

    def panstarrs_query (self, ra_deg, dec_deg, ra_max_deg,dec_max_deg, mindet=5, maxsources=5000):

        returnTable = None
        input = None


        if self.indexTable is not None:
            # Fan out the search into segments as defined in indexTable


            needed_catalogs =  (self.indexTable['R_MIN'] <= ra_max_deg)\
                              & (self.indexTable['R_MAX'] > ra_deg)\
                              & (self.indexTable['D_MIN'] <= dec_max_deg)\
                              & (self.indexTable['D_MAX'] > dec_deg)

            #print (needed_catalogs)

            first = True
            f= open ("tmp.cat", "w+")
            for needed in self.indexTable[needed_catalogs]:
                #print (needed)
                cachename = os.path.join(self.localCacheIndexBaseDir, needed['NAME'])
                print ("Reqeusted area falls into segment: %s" % cachename)


                iteminput = self.do_panstarrs_query(needed['R_MIN'], \
                                                    needed['D_MIN'], \
                                                    needed['R_MAX'], \
                                                    needed['D_MAX'], \
                                                    mindet=mindet, maxsources=5000, cacheFile = cachename)



                if iteminput is not None:
                    if not first:
                        iteminput.readline()
                        iteminput.readline()
                    f.write(iteminput.read(None))
                    iteminput.close()
                    first = False
            f.close()
            input = "tmp.cat"
        else:
            # no Index defined, go ahead and do a streight query to PS.
            # TODO: check for 360 degree roll-over in RA
            input = self.do_panstarrs_query(ra_deg, dec_deg, ra_max_deg,dec_max_deg, mindet=mindet,
                    maxsources=maxsources)

        if input is not None:



            try:
                returnTable = np.genfromtxt (input, names=True, skip_header=1, comments='#', delimiter=',', invalid_raise=False)
                print ("Got  lines from query: " , returnTable.shape)
            except:
                returnTable = None

            self.PS1toSDSS(returnTable)

        return returnTable


    def do_panstarrs_query(self, ra_deg, dec_deg, ra_max_deg,dec_max_deg, mindet=5,
                    maxsources=500000,
                    server=('http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx'),
                    cacheFile = None):

        """ Online query to PS1 from MAST with check if data cache exists at given location

            Returns an object taht is radable by np.genfromtxt
        """

        #print ((cacheFile is not None) , (os.path.isfile(cacheFile)))

        if (cacheFile is not None) and (os.path.isfile(cacheFile)):
            # if lcoal cache exists, use it.
            # Sounds trivial, but: a chacehFile should contain the ascii text and will be readilyl readable
            # numpy. Nothing to do here.
            print ("Detected existing valid cache file. Using it")
            with open (cacheFile, 'r') as f:
                ret= StringIO.StringIO(f.read())
                f.close()
                return ret
            return StringIO.StringIO()
        else:
            # No cache provided, need to retrieve fresh data. Upon retrieval, data to be stored in cache
            area = ("%6f,%6f,%6f,%6f" % (ra_deg,dec_deg,ra_max_deg,dec_max_deg))
            print (area)
            print ("Start retrieve catalog data at ra % 8.4f dec % 8.4f % 8.4f %8.4f" % (ra_deg,dec_deg, ra_max_deg, dec_max_deg))

            try:
             r = requests.get(server,
                    params= {'CAT':"PS1V3OBJECTS", 'BBOX':area , 'MAXOBJ': maxsources,
                 'FORMAT': 'CSV',
                    'ndetections': ('>%d' % mindet)})
            except:
                r=None

            if (r is None):
                return


            # need at least a few lines of input to reasonably proceed.

            nlines = r.text.count('\n')
            if (nlines > 10):
                # Let's write the local output into a cache file.

                if (cacheFile != None):

                    if not os.path.isfile(cacheFile):
                        try:
                            # create parent directory if needed.
                            if not os.path.exists(os.path.dirname(cacheFile)):
                                os.makedirs(os.path.dirname(cacheFile))

                            # write out file.
                            with open(cacheFile, "w") as f:
                                print ("Writing local PS1 cache file: %s" % (cacheFile))
                                f.write(r.text)
                                f.close()

                        except Exception as error: # Guard against race condition
                            print ("Failure while creating cache file directory: " + error.message)


                    else:
                        print ("Cached file already exists, not writing again")




                input = StringIO.StringIO(r.text)
                return input
            else:
                print ("Did not retireve enough lines from query!\n %s \n\n%s" % (r.url,r.text))


        print ("\t failed retrieve table)")
        return None



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

# Example query

from astropy import units as u

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import re



import glob

inputlist = glob.glob("../testing/g/*.fits.fz")


photzpStage = PhotCalib()


#ps1 = PS1Catalog(None)

for image in inputlist:
    print ("Work in image %s" % image)
    photzpStage.analyzeImage(image)
