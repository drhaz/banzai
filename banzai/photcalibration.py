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
import sys
import logging

_logger = logging.getLogger(__name__)

__author__ = 'dharbeck'


class PhotCalib():
    # To be replaced with map:
    # LCO filter -> sdss filter name, sdsss g-i color term, airmass term, default zero point.
    # possibly need this for all site:telescope:camera:filter settings. Maybe start with one
    # default and see where it goes.
    ps1 = None

    def __init__(self, ps1dir):
        # super(PhotCalib, self).__init__(pipeline_context)
        # TODO: capture error if no ps1 catalog is found at location
        self.ps1 = PS1IPP(ps1dir)

    def isInCatalogFootprint(self, ra, dec):
        """ Verify if image is in catalog footprint.
            TODO: Account for image field of view
        """

        # PanSTARRS has valid entries for DEc > - 30 degrees
        return dec >= -30.0

    def do_stage(self, images):
        """ future hook for BANZAI pipeline integration

        """

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

        if not self.isInCatalogFootprint(ra, dec):
            _logger.debug("Image not in PS1 footprint. Ignoring")
            testimage.close()
            return None

        retCatalog['exptime'] = testimage['SCI'].header['EXPTIME']

        # Safeguard against a division by zero downstream.
        if (retCatalog['exptime'] <= 0):
            _logger.debug("Image has exposure time <= 0 % 8.2f. Ignoring." % (retCatalog['exptime']))
            testimage.close()
            return None

        retCatalog['instfilter'] = testimage['SCI'].header['FILTER']
        retCatalog['airmass'] = testimage['SCI'].header['AIRMASS']
        retCatalog['dateobs'] = testimage['SCI'].header['DATE-OBS']
        retCatalog['instrument'] = testimage['SCI'].header['INSTRUME']
        retCatalog['siteid'] = testimage['SCI'].header['SITEID']
        retCatalog['domid'] = testimage['SCI'].header['ENCLOSUR']
        retCatalog['telescope'] = testimage['SCI'].header['TELID']
        retCatalog['FOCOBOFF'] = testimage['SCI'].header['FOCOBOFF']

        # Check if filter is supported
        if retCatalog['instfilter'] not in self.ps1.FILTERMAPPING:
            _logger.debug("Filter %s not viable for photometric calibration. Sorry" % (retCatalog['instfilter']))
            testimage.close()
            return None

        # Check if exposure time is long enough
        if (retCatalog['exptime'] < 60):
            _logger.debug("Exposure %s time is deemed too short, ignoring" % (retCatalog['exptime']))
            testimage.close()
            return None

        # verify there is no deliberate defocus
        if (retCatalog['FOCOBOFF'] is not None) and (retCatalog['FOCOBOFF'] != 0):
            _logger.debug("Exposure is deliberately defocussed by %s, ignoring" % (retCatalog['FOCOBOFF']))
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
            testimage.close()
            return None

        # Transform the image catalog to RA / Dec based on the WCS solution in the header.
        # TODO: rerun astrometry.net with a higher order distortion model

        image_wcs = WCS(testimage['SCI'].header)
        try:
            ras, decs = image_wcs.all_pix2world(instCatalog['x'], instCatalog['y'], 1)
        except:
            _logger.error("Failed to convert images coordinates to world coordinates. Giving up on file.")
            testimage.close()
            return None

        # Now we have all we wanted from the input image, close it
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

    def reject_outliers(self, data, m=2):
        """
        Reject data from vector that are > m std deviations from median
        :param m:
        :return:
        """
        std = np.std (data)
        return data[abs(data - np.median(data)) < m * std]

    def analyzeImage(self, imageName, pickle="photzp.db",
                     outputimageRootDir=None, sqlite3cur = None):
        """ Do full photometric zeropoint analysis on an image
         """

        # TODO: Make this procedure thread safe so it can be accelerated a bit.

        retCatalog = self.loadFitsCatalog(imageName)

        if (retCatalog is None) or (retCatalog['instmag'] is None) or (len(retCatalog['ra']) < 10):
            _logger.debug("Not enough stars to fit %s" % (imageName,))
            return

        # calculate the per star zeropoint
        magZP = retCatalog['refmag'] - retCatalog['instmag']

        refmag = retCatalog['refmag']
        ra = retCatalog['ra']
        dec = retCatalog['dec']
        refcol = retCatalog['refcol']

        # Calculate the photometric zeropoint.
        # TODO: Robust median w/ rejection, error propagation.

        cleandata =   self.reject_outliers(magZP, 3)
        photzp = np.median(cleandata)
        photzpsig = np.std (cleandata)

        # calculate color term

        try:
            cond = (refcol > 0) & (refcol < 3) & (np.abs(magZP - photzp) < 0.75)
            colorparams = np.polyfit(refcol[cond], (magZP - photzp)[cond], 1)
            color_p = np.poly1d(colorparams)
            delta = np.abs(magZP - photzp - color_p(refcol))
            cond = (delta < 0.2)
            colorparams = np.polyfit(refcol[cond], (magZP - photzp)[cond], 1)
            color_p = np.poly1d(colorparams)
            colorterm = colorparams[0]

        except:
            _logger.warn("could not fit a color term. ")
            color_p = None
            colorterm = 0

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
            plt.savefig("%s/%s_%s_zp.png" % (outputimageRootDir, outbasename, retCatalog['instfilter']))
            plt.close()

            plt.figure()
            plt.plot(refcol, magZP - photzp, '.')
            if color_p is not None:
                xp = np.linspace(-0.5, 3.5, 10)
                plt.plot(xp, color_p(xp), '-', label="color term fit % 6.4f" % (colorterm))
                plt.legend()
            # plt.xlim([-0.5, 3])
            # plt.ylim([-1, 1])
            plt.xlabel("(g-r)_{SDSS} Reference")
            plt.ylabel("Reference Mag - Instrumnetal Mag - ZP (%5.2f) %s" % (photzp, retCatalog['instfilter']))
            plt.title("Color correction %s " % (outbasename))
            plt.savefig("%s/%s_%s_color.png" % (outputimageRootDir, outbasename, retCatalog['instfilter']))
            plt.close()

            # plt.figure()
            # plt.scatter(ra, dec, c=magZP - photzp, vmin=-0.2, vmax=0.2, edgecolor='none',
            #             s=9, cmap=matplotlib.pyplot.cm.get_cmap('nipy_spectral'))
            # plt.colorbar()
            # plt.title("Spacial variation of phot. Zeropoint %s %s" % (outbasename,retCatalog['instfilter']))
            # plt.xlabel("RA")
            # plt.ylabel("Dec")
            # plt.savefig("%s/%s_%s_zpmap.png" % (outputimageRootDir, outbasename,retCatalog['instfilter']))
            # plt.close()

        # TODO: Make this thread safe, e.g., write to transactional database, or return values for storing externally.
        with open(pickle, 'a') as f:
            output = "%s %s %s %s %s %s %s %s % 6.3f  % 6.3f  % 6.3f\n" % (
                imageName, retCatalog['dateobs'], retCatalog['siteid'], retCatalog['domid'],
                retCatalog['telescope'], retCatalog['instrument'], retCatalog['instfilter'],
                retCatalog['airmass'], photzp, colorterm, photzpsig)
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
    ps1colorterms['g'] = [-0.01808, -0.13595, +0.01941, -0.00183][::-1]
    ps1colorterms['r'] = [-0.01836, -0.03577, +0.02612, -0.00558][::-1]
    ps1colorterms['i'] = [+0.01170, -0.00400, +0.00066, -0.00058][::-1]
    ps1colorterms['z'] = [-0.01062, +0.07529, -0.03592, +0.00890][::-1]

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

        # A lot of safeguarding boiler plate to ensure catalog files are valid.
        if (self.basedir is None) or (not os.path.isdir(self.basedir)):
            _logger.error("Unable to find reference catalog: %s" % (str(self.basedir)))
            return None

        # Load the SkyTable so we know in what files to look for the catalog"
        _logger.debug("Using catalog found in %s" % (self.basedir))
        skytable_filename = "%s/SkyTable.fits" % (self.basedir)
        if (not os.path.isfile(skytable_filename)):
            _logger.error("Unable to find catalog index file in %s!" % (self.basedir))
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

        _logger.debug("Querying catalog: Ra=%f...%f Dec=%f...%f" % (min_ra, max_ra, min_dec, max_dec))

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

        _logger.debug("# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec))

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
        _logger.debug(files_to_read)

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
                    _logger.warning(
                        "Catalog file (%s) not found (base-dir: %s)" % (os.path.abspath(catalogfile), self.basedir))
                    continue

            try:
                hdu_cat = fits.open(catalogfile)
            except:
                _logger.warning("Unable to open catalog file %s" % (catalogfile))
                continue

            catalog_filenames.append(catalogfile)
            _logger.debug("Adding %s to list of catalog files being used" % (catalogfile))

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
            _logger.debug("Read %d sources from %s" % (array_to_add.shape[0], catalogname))

            if (full_catalog is None):
                full_catalog = array_to_add
            else:
                full_catalog = np.append(full_catalog, array_to_add, axis=0)
                # print photom_grizy[:3,:]

            if (full_catalog is None):
                _logger.warning("No stars found in area %s, %s from catalog %s" % (
                    str(ra), str(dec),
                    # ra[0], ra[1], dec[0], dec[1],
                    self.basedir))
            else:
                _logger.debug(
                    "Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read)))

        self.PS1toSDSS(full_catalog)
        return full_catalog


def crawlDirectory (directory, imagedb,  args):


    search = "%s/*-[es]91.fits.fz" % (directory)
    inputlist = glob.glob(search)
    _logger.debug("Found %d entries. Cleaning duplicate entries..." % len(inputlist))
    # TODO: Do not lie, and actually do clean for duplicate entries!
    _logger.debug("Starting analysis")


    photzpStage = PhotCalib(args.ps1dir)
    for image in inputlist:
        image = image.rstrip()
        photzpStage.analyzeImage(image, pickle=imagedb, outputimageRootDir=args.outputimageRootDir)

def crawlSiteCameraArchive(site, camera, args, date=None):
    '''
    Process in the archive

    :param site:
    :param camera:
    :param args:
    :param date:
    :return:
    '''

    if date is None:
        date = '*'

    if site is None:
        site = '*'

    imagedb = "%s/%s-%s.db" % (args.imagedbPrefix, site, camera)
    searchdir = "%s/%s/%s/%s/processed" % (args.rootdir, site, camera, date)
    #search = "%s/%s/%s/%s/preview/*-[es]11.fits.fz" % (args.rootdir, site, camera, date)
    _logger.info("File search string is: %s" % (searchdir))

    crawlDirectory(searchdir, imagedb, args)




def crawlSite(site, type, args):
    """ Search for all cameras of a given type (fl, kb,fs) in a site directory and processed them """
    searchdir = "%s/%s/%s*" % (args.rootdir, site, type)

    cameralist = glob.glob(searchdir)
    cameras = []
    for candidate in cameralist:
        cameras.append((site, os.path.basename(os.path.normpath(candidate))))

    for setup in cameras:
        print(setup[0], setup[1])
        crawlSiteCameraArchive(setup[0], setup[1], args, args.date)


def parseCommandLine():
    """ Read command line parameters
    """

    parser = argparse.ArgumentParser(
        description='Determine photometric zeropoint of banzai-reduced LCO imaging data.')

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO'],
                        help='Set the debug level')

    parser.add_argument('--ps1dir', dest='ps1dir', default='~/Catalogs/ps1odi/panstarrs/',
                        help='Directory of PS1 catalog')

    parser.add_argument("--diagnosticplotsdir", dest='outputimageRootDir', default=None,
                        help='Output directory for diagnostic photometry plots. No plots generated if option is omitted. ')

    parser.add_argument('--imagedbPrefix', dest='imagedbPrefix', default='~/lcozpplots',
                        help='Result output directory. .db file is written here')

    parser.add_argument('--imagerootdir', dest='rootdir', default='/nfs/archive/engineering',
                        help="LCO archive root directory")

    parser.add_argument('--site', dest='site', default=None, help='sites code for camera')

    parser.add_argument('--date', dest='date', default=None, help='Specific date to process.')

    cameragroup = parser.add_mutually_exclusive_group()
    cameragroup.add_argument('--camera', dest='camera', default=None, help='specific camera to process. ')
    cameragroup.add_argument('--cameratype', dest='cameratype', default=None, choices=['fs', 'fl', 'kb'],
                             help='camera type to process at selected sites to process. ')
    cameragroup.add_argument('--crawldirectory', default=None, type=str, help="process all reduced image in specific directoy")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    args.imagedbPrefix = os.path.expanduser(args.imagedbPrefix)
    if args.outputimageRootDir is not None:
        args.outputimageRootDir = os.path.expanduser(args.outputimageRootDir)

    if args.crawldirectory is not None:
        args.crawldirectory = os.path.expanduser(args.crawldirectory)

    args.ps1dir = os.path.expanduser(args.ps1dir)

    return args


if __name__ == '__main__':

    args = parseCommandLine()

    if args.cameratype is not None:

        cameras = [camera for camera in args.cameratype.split(',')]

        if args.site is not None:
            sites = [site for site in args.site.split(',')]
        else:
            sites = ('lsc', 'cpt', 'ogg', 'coj', 'tfn', 'elp')

        print("Crawling through camera types ", cameras, " at sites ", sites, " for date ", args.date)
        for site in sites:
            for cameratype in cameras:
                crawlSite(site, cameratype, args)

    elif args.camera is not None:

        if args.site is None:
            sites = '*'
        else:
            sites = args.site

        print("Calibrating camera ", args.camera, " at site ", sites, ' for date ', args.date)
        crawlSiteCameraArchive(sites, args.camera, args, date=args.date)


    elif args.crawldirectory is not None:
        imagedb = "%s/%s" % (args.crawldirectory, 'imagezp.db')

        crawlDirectory(args.crawldirectory, imagedb, args)


    else:
        print("Need to specify either a camera, or a camera type.")
    sys.exit(0)
