from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import argparse
import re
import glob
import os
import math

import logging

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)
logging.basicConfig(format='%(asctime)-15s %(levelname)s\t %(message)s', level=logging.NOTSET)


__author__ = 'dharbeck'


class PhotCalib():
    # To be replaced with map:
    # LCO filter -> sdss filter name, sdsss g-i color term, airmass term, default zero point.
    # possibly need this for all site:telescope:camera:filter settings. Maybe start with one
    # default and see where it goes.
    ps1 = None

    def __init__(self):
        # super(PhotCalib, self).__init__(pipeline_context)
        self.ps1 = PS1IPP("/home/dharbeck/Catalogs/ps1odi/panstarrs/")

    def isInCatalogFootprint(dec):
        """ Verify if image is in catalog footprint.
            TODO: Account for image field of view
        """

        # PanSTARRS has valid entries for DEc > - 30 degrees
        return dec >= -30.0

    def do_stage(self, images):

        for i, image in enumerate(images):
            pass
            # logging_tags = logs.image_config_to_tags(image, self.group_by_keywords)

    def loadFitsCatalog(self, image):
        """ Load the banzai-generated photometry catalog from  'CAT' extension, queries PS1 catalog for image FoV, and
        returns a cross-matched catalog.

        Errors conditions:
         if photometricfilter is not supported, none is returned
         if exposure time is < 60 seconds, None is returned.



        :param image: input fits image path+name
        :return:
        """
        # Build up a baseline refernce catlog with meta data
        retCatalog = {'fname': image,
                      'instmag': None
                      }

        # Read banzai star catalog
        testimage = fits.open(image)

        # Boilerplate grab of status infomration
        ra = testimage['SCI'].header['CRVAL1']
        dec = testimage['SCI'].header['CRVAL2']
        retCatalog['exptime'] = testimage['SCI'].header['EXPTIME']

        # Safeguard against a division by zero downstream.
        if (retCatalog['exptime'] <= 0):
            retCatalog['exptime'] = 1

        retCatalog['instfilter'] = testimage['SCI'].header['FILTER']
        retCatalog['airmass'] = testimage['SCI'].header['AIRMASS']
        retCatalog['dateobs'] = testimage['SCI'].header['DATE-OBS']
        retCatalog['instrument'] = testimage['SCI'].header['INSTRUME']
        retCatalog['siteid'] = testimage['SCI'].header['SITEID']
        retCatalog['domid'] = testimage['SCI'].header['ENCLOSUR']
        retCatalog['telescope'] = testimage['SCI'].header['TELID']

        # Check if filter is supported
        if retCatalog['instfilter'] not in self.ps1.FILTERMAPPING:
            _logger.debug("Filter not viable for photometric calibration. Sorry")
            testimage.close()
            return None

        if (retCatalog['exptime'] < 60):
            _logger.debug("Exposure %s time is deemed too short, ignoring" % (retCatalog['exptime']))
            testimage.close()
            return None

        # Get the instrumental filter and the matching reference catalog filter names.
        referenceInformation = self.ps1.FILTERMAPPING[retCatalog['instfilter']]
        referenceFilterName = referenceInformation['refMag']

        # Load photometry catalog from image, and transform into RA/Dec coordinates
        try:
            instCatalog = testimage['CAT'].data
        except:
            _logger.error("No extension \'CAT\' available for image %s, skipping." % (image))
            return None

        # Transform the image catalog to RA / Dec based on the WCS solution in the header.
        # TODO: rerun astrometry.net with a higher order distortion model

        image_wcs = WCS(testimage['SCI'].header)
        try:
            ras, decs = image_wcs.all_pix2world(instCatalog['x'], instCatalog['y'], 1)
        except:
            _logger.error("Failed to convert images ccordinates to world coordinates. Giving up on file.")
            return None

        # Now we have all we waned from the input image, close it
        testimage.close()

        # Query reference catalog
        reftable = self.ps1.get_reference_catalog(ra, dec, 0.25)
        if reftable is None:
            _logger.warn("Failure on image %s, no reference catalog received." % (image))
            return None

        # Start the catalog matching, using astropy skycoords built-in functions.
        cInstrument = SkyCoord(ra=ras * u.degree, dec=decs * u.degree)
        cReference = SkyCoord(ra=reftable['RA'] * u.degree, dec=reftable['DEC'] * u.degree)
        idx, d2d, d3d = cReference.match_to_catalog_sky(cInstrument)

        # Reshuffle the source catalog to index-match the reference catalog.
        # There is probably as smarter way of doing this!
        instCatalogT = np.transpose(instCatalog)[idx]
        instCatalog = np.transpose(instCatalogT)

        # Measure the distance between matched pairs. Important to down-select viable pairs.
        distance = cReference.separation(cInstrument[idx]).arcsecond

        # Define a reasonable condition on what is a good match on good photometry
        condition = (distance < 5) & (instCatalog['FLUX'] > 0) & (reftable[referenceFilterName] > 0) & (
            reftable[referenceFilterName] < 26)

        # Calculate instrumental magnitude from PSF instrument photometry
        instmag = -2.5 * np.log10(instCatalog['FLUX'][condition] / retCatalog['exptime'])

        # Calculate the magnitude difference between reference and inst catalog
        retCatalog['instmag'] = instmag
        retCatalog['refcol'] = (reftable['g'] - reftable['i'])[condition]

        retCatalog['refmag'] = reftable[referenceFilterName][condition]
        retCatalog['ra'] = reftable['RA'][condition]
        retCatalog['dec'] = reftable['DEC'][condition]
        retCatalog['matchDistance'] = distance[condition]
        # TODO: Read photometric error columns from reference and instrument catalogs, properly propagate error.

        return retCatalog

    def analyzeImage(self, imageName, pickle="photzp.db", outputimageRootDir=None):
        """ Do full photometric zeropoint analysis on an image"""

        retCatalog = self.loadFitsCatalog(imageName)

        if (retCatalog is None) or (retCatalog['instmag'] is None):
            return

        # calculate the per star zeropoint
        magZP = retCatalog['refmag'] - retCatalog['instmag']

        refmag = retCatalog['refmag']
        ra = retCatalog['ra']
        dec = retCatalog['dec']
        refcol = retCatalog['refcol']

        # Calculate the photometric zeropoint.
        # TODO: Robust median w/ rejection, error propagation.
        photzp = np.median(magZP)

        if (outputimageRootDir is not None) and (os.path.exists(outputimageRootDir)):
            outbasename = os.path.basename(imageName)
            outbasename = re.sub('.fits.fz', '', outbasename)
            plt.figure()
            plt.plot(refmag, magZP, '.')
            plt.xlim([10, 22])
            plt.ylim([20, 26])
            plt.axhline(y=photzp, color='r', linestyle='-')
            plt.xlabel("Reference catalog mag")
            plt.ylabel("Reference Mag - Instrumnetal Mag (%s)" % (retCatalog['instfilter']))
            plt.title("Photometric zeropoint %s %5.2f" % (outbasename, photzp))
            plt.savefig("%s/%s_%s_zp.png" % (outputimageRootDir,outbasename,retCatalog['instfilter']))
            plt.close()

            plt.figure()
            plt.plot(refcol, magZP - photzp, '.')
            plt.xlim([-0.5, 3])
            plt.ylim([-1, 1])
            plt.xlabel("(g-r)_{SDSS} Reference")
            plt.ylabel("Reference Mag - Instrumnetal Mag - ZP (%5.2f) %s" % (photzp,retCatalog['instfilter']))
            plt.title("Color correction %s " % (outbasename))
            plt.savefig("%s/%s_%s_color.png" % (outputimageRootDir,outbasename,retCatalog['instfilter']))
            plt.close()

            plt.figure()
            plt.scatter(ra, dec, c=magZP - photzp, vmin=-0.2, vmax=0.2, edgecolor='none',
                        s=9, cmap=matplotlib.pyplot.cm.get_cmap('nipy_spectral'))
            plt.colorbar()
            plt.title("Spacial variation of phot. Zeropoint %s %s" % (outbasename,retCatalog['instfilter']))
            plt.xlabel("RA")
            plt.ylabel("Dec")
            plt.savefig("%s/%s_%s_zpmap.png" % (outputimageRootDir, outbasename,retCatalog['instfilter']))
            plt.close()


        with open(pickle, 'a') as f:
            output = "%s %s %s %s %s %s %s %s % 6.3f \n" % (
                imageName, retCatalog['dateobs'], retCatalog['siteid'], retCatalog['domid'],
                retCatalog['telescope'], retCatalog['instrument'], retCatalog['instfilter'],
                retCatalog['airmass'], photzp)
            _logger.info(output)
            f.write(output)
            f.close()

        return photzp




class PS1IPP:
    """ Class to access local, distilled copy of PS1 data release.

        Based on code from WIYN ODI quickreduce pipeline, developed by Ralf Kotula. See:
        https://github.com/WIYN-ODI/QuickReduce
    """

    FILTERMAPPING = {}
    FILTERMAPPING['gp'] = {'refMag': 'g', 'colorTerm': 0.0, 'airmassTerm': 0.20, 'defaultZP': 0.0}
    FILTERMAPPING['rp'] = {'refMag': 'r', 'colorTerm': 0.0, 'airmassTerm': 0.12, 'defaultZP': 0.0}
    FILTERMAPPING['ip'] = {'refMag': 'i', 'colorTerm': 0.0, 'airmassTerm': 0.08, 'defaultZP': 0.0}
    FILTERMAPPING['zp'] = {'refMag': 'z', 'colorTerm': 0.0, 'airmassTerm': 0.05, 'defaultZP': 0.0}

    ###  PS to SDSS color transformations according to  Finkbeiner 2016
    ###  http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/meta#apj522061s2-4 Table 2
    ###  Note that this transformation is valid for stars only. For the purpose of photometric
    ###  calibration, it is desirable to select point sources only from the input catalog.

    ## Why reverse the order of the color term entries? Data are entered in the order as they are
    ## shown in paper. Reverse after the fact to avoid confusion when looking at paper

    ps1colorterms = {}
    ps1colorterms['g'] = [-0.01808, -0.13595, 0.01941, -0.00183][::-1]
    ps1colorterms['r'] = [-0.01836, -0.03577, 0.02612, -0.00558][::-1]
    ps1colorterms['i'] = [0.01170, -0.00400, 0.00066, -0.00058][::-1]
    ps1colorterms['z'] = [-0.01062, 0.07529, -0.03592, 0.00890][::-1]

    def __init__(self, basedir):
        self.basedir = basedir
        self.skytable = None

    def PS1toSDSS(self, table):
        """
        Modify table in situ from PS1 to SDSS, requires column names compatible with ps1colorterms definition.

        :param table:
        :return: modified table.
        """
        if table is not None:
            pscolor = table['g'] - table['i']
            for filter in self.ps1colorterms:
                colorcorrection = np.polyval(self.ps1colorterms[filter], pscolor)
                table[filter] -= colorcorrection

        return table

    def get_reference_catalog(self, ra, dec, radius, overwrite_select=False):
        """ Read i fits table from local catalog copy. Concatenate tables columns
           from different fits tables for full coverage.
        """

        # TODO: integrate into logging schema
        logger = logging

        # A lot of safeguarding boiler plate to ensure catalog files are valid.
        if (self.basedir is None) or (not os.path.isdir(self.basedir)):
            logger.warning("Unable to find reference catalog: %s" % (str(self.basedir)))
            return None

        # Load the SkyTable so we know in what files to look for the catalog"
        logger.debug("Using catalog found in %s" % (self.basedir))
        skytable_filename = "%s/SkyTable.fits" % (self.basedir)
        if (not os.path.isfile(skytable_filename)):
            logger.error("Unable to find catalog index file in %s!" % (self.basedir))
            return None

        # Read in the master index hdu
        skytable_hdu = fits.open(skytable_filename)
        skytable = skytable_hdu['SKY_REGION'].data

        # Select entries that match our list
        # print ra, dec, radius, type(ra), type(dec), type(radius)
        # logger.debug("# Searching for stars within %.1f degress around %f , %f ..." % (radius, ra, dec))

        if (not radius == None and radius > 0):
            min_dec = dec - radius
            max_dec = dec + radius
            min_ra = ra - radius / math.cos(math.radians(dec))
            max_ra = ra + radius / math.cos(math.radians(dec))
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

        logger.debug("# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec))

        try:
            needed_catalogs = (skytable['PARENT'] > 0) & (skytable['PARENT'] < 25) & \
                              (skytable['R_MAX'] > min_ra) & (skytable['R_MIN'] < max_ra) & \
                              (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)
        except KeyError:
            # try without the PARENT field
            needed_catalogs = (skytable['R_MAX'] > min_ra) & (skytable['R_MIN'] < max_ra) & \
                              (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)

        # print skytable[needed_catalogs]

        files_to_read = skytable['NAME'][needed_catalogs]
        files_to_read = [f.strip() for f in files_to_read]
        logger.debug(files_to_read)

        skytable_hdu.close()  # Warning: might erase the loaded data, might need to copy array!

        # Now quickly go over the list and take care of all filenames that still have a 0x00 in them
        for i in range(len(files_to_read)):
            found_at = files_to_read[i].find('\0')
            if (found_at > 0):
                files_to_read[i] = files_to_read[i][:found_at]

        # Load all frames, one by one, and select all stars in the valid range.
        # Then add them to the catalog with RAs and DECs
        full_catalog = None  # numpy.zeros(shape=(0,6))
        catalog_filenames = []

        # Start iterating though catalogs
        for catalogname in files_to_read:

            catalogfile = "%s/%s" % (self.basedir, catalogname)
            # print catalogfile
            if (not os.path.isfile(catalogfile)):
                # not a file, try adding .fits to the end of the filename
                if (os.path.isfile(catalogfile + ".fits")):
                    catalogfile += ".fits"
                else:
                    # neither option (w/ or w/o .fits added is a file)
                    logger.warning(
                        "Catalog file (%s) not found (base-dir: %s)" % (os.path.abspath(catalogfile), self.basedir))
                    continue

            try:
                hdu_cat = fits.open(catalogfile)
            except:
                logger.warning("Unable to open catalog file %s" % (catalogfile))
                continue

            catalog_filenames.append(catalogfile)
            logger.debug("Adding %s to list of catalog files being used" % (catalogfile))

            # read table into a nd-array buffer
            cat_full = hdu_cat[1].data
            hdu_cat.close()

            # Read the RA and DEC values
            cat_ra = cat_full['RA']
            cat_dec = cat_full['DEC']

            # To select the right region, shift a temporary catalog
            cat_ra_shifted = cat_ra
            if (max_ra > 360.):
                cat_ra_shifted[cat_ra < 180] += 360
            elif (min_ra < 0):
                cat_ra_shifted[cat_ra > 180] -= 360

            select_from_cat = (cat_ra_shifted > min_ra) & (cat_ra_shifted < max_ra) & (cat_dec > min_dec) & (
                cat_dec < max_dec)

            array_to_add = cat_full[select_from_cat]
            logger.debug("Read %d sources from %s" % (array_to_add.shape[0], catalogname))

            if (full_catalog is None):
                full_catalog = array_to_add
            else:
                full_catalog = np.append(full_catalog, array_to_add, axis=0)
                # print photom_grizy[:3,:]

            if (full_catalog is None):
                logger.warning("No stars found in area %s, %s from catalog %s" % (
                    str(ra), str(dec),
                    # ra[0], ra[1], dec[0], dec[1],
                    self.basedir))
            else:
                logger.debug("Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read)))

        self.PS1toSDSS(full_catalog)
        return full_catalog

#import longtermphotzp

def crawlDirectory (site, camera, basedir='/nfs/archive/engineering',outputimageRootDir=None):

    photzpStage = PhotCalib()

    imagedb = "/home/dharbeck/lcozpplots/%s-%s.db" % (site,camera)
    search = "%s/%s/%s/*/processed/*-e91.fits.fz" % (basedir, site, camera)
    inputlist = glob.glob(search)

    print("Found %d entries. Cleaning duplicate entries..." % len(inputlist))
    print("Starting analysis")

    for image in inputlist:
        image = image.rstrip()
        photzpStage.analyzeImage(image, pickle=imagedb, outputimageRootDir=outputimageRootDir)


def crawlSite ( site, type, basedir='/nfs/archive/engineering', outputimageRootDir=None):
    searchdir = "%s/%s/%s*" % (basedir, site, type)

    cameralist = glob.glob (searchdir)
    cameras = []
    for candidate in cameralist:
        cameras.append ( (site, os.path.basename(os.path.normpath(candidate)), outputimageRootDir))


    for setup in cameras:
        print ( setup[0],setup[1],setup[2])
        #crawlDirectory(setup[0],setup[1],setup[2])



def parseCommandLine ():
    parser = argparse.ArgumentParser(
        description='Determine photometric zeropoint of banzai-reduced LCO imaging data.')
    parser.add_argument("--diagnosticplotsdir", dest='outputimageRootDir', default=None,
                    help='Output directory for diagnostic plots. No plots generated if option is omitted. ')

    parser.add_argument('--imagedbPrefix', dest='imagedbPrefix', default='/home/dharbeck/lcozpplots', help = 'Result output directory')
    parser.add_argument('--rootdir', dest='rootdir', default='/nfs/archive/engineering', help="data root directory")
    parser.add_argument('--site', dest='site', default='lsc', help='sites code for camera')
    parser.add_argument('--camera', dest='camera', default='fl03', help='camera to process')


    args = parser.parse_args()
    print (args)
    pass

import multiprocessing
if __name__ == '__main__':
   # crawlDirectory('elp','fl05')
    global args

    args = parseCommandLine()

    sites = ('lsc','cpt','ogg','coj','tfn')
    cameras = ("fl","fs","kb")



    for site in sites:
        for cameratype in cameras:
            crawlSite(site, cameratype)

