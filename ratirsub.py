import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator, time
import pyfits
from types import *
import numpy as np

from iqutils import *
from iqpkg import *

# Necessary packages
iraf.images()
iraf.immatch()
iraf.imfilter()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.digiphot()
iraf.apphot()

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets
imcombine=iraf.imcombine

pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

globclob=yes
globver=yes

###########################################################################

def p60hotpants(inlis, refimage, outimage, tu=50000, iu=50000,
                stamps=None, nsx=4, nsy=4, ko=0, bgo=0, radius=10, 
                tlow=0, ilow=0, sthresh=5.0, ng=None,
                ngref=[3, 6, 0.70, 4, 1.50, 2, 3.00], scimage=no,
                rms=True):

    '''P60 Subtraction using HOTPANTS'''

    subimages=iraffiles(inlis)

    for image in subimages:

        root = image.split('.')[0]
        scmd="hotpants -inim %s -tmplim %s -outim %s -tu %.2f -tuk %.2f -iu %.2f -iuk %.2f -savexy %s.st -ko %i -bgo %i -nsx %i -nsy %i -r %i -rss %i -tl %f -il %f -ft %f -v 0 -c t -n i" % (image, refimage, outimage, tu, tu, iu, iu, root, ko, bgo, nsx, nsy, radius, radius*1.5, tlow, ilow, sthresh)
        if not (stamps==None):
            scmd += " -ssf %s -afssc 0" % stamps 
        if (ng==None):
            scmd += " -ng %i " % ngref[0]
            seepix=get_head(image,'SEEPIX')
            for k in range(1, len(ngref)-1, 2):
                scmd+="%i %.2f " % (ngref[k], ngref[k+1]*float(seepix)/3.0)
        if (scimage==yes):
            scmd += ' -n i '
        if (rms==True) and os.path.exists("%s.rms.fits" % refimage[:-5]):
        	scmd += ' -tni %s.rms.fits ' % refimage[:-5]
        if (rms==True) and os.path.exists("%s.rms.fits" % image[:-5]):
        	scmd += ' -ini %s.rms.fits ' % image[:-5]
        if os.path.exists("%s.mask.fits" % refimage[:-5]):
            scmd += ' -tmi %s.mask.fits ' % refimage[:-5]
        if os.path.exists("%s.mask.fits" % image[:-5]):
            scmd += ' -imi %s.mask.fits ' % image[:-5] 
        print scmd
        cmd=os.popen(scmd,'r')
        hlines=cmd.readlines()
        cmd.close()

    return

############################################################################

def p60swarp(image, outfile, ractr=None, dcctr=None, pixscale=None, size=None,
             backsub=no, refimage=None, rms=True, mask=True):

    '''Run SWarp on P60 images'''

    swarpcmd='swarp %s ' % image
    if ractr!=None or dcctr!=None:
        swarpcmd+='-CENTER_TYPE MANUAL -CENTER "%s, %s" ' % (ractr, dcctr)
    if pixscale!=None:
        swarpcmd+='-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %f ' % pixscale
    if size!=None:
        swarpcmd+='-IMAGE_SIZE "%i, %i" ' % (size[0], size[1])
    if backsub==no:
        swarpcmd+='-SUBTRACT_BACK N '
    swarpcmd+='-COPY_KEYWORDS OBJECT,SKYSUB,SKYBKG,SEEPIX '

    if refimage!=None:
		[nax1,nax2,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2] = get_head(refimage, ["NAXIS1", "NAXIS2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"])
		hout = open("coadd.head", "w")
		hout.write("CRVAL1  =     %.8g\n" % crval1)
		hout.write("CRVAL2  =     %.8g\n" % crval2)
		hout.write("CRPIX1  =     %.8g\n" % crpix1)
		hout.write("CRPIX2  =     %.8g\n" % crpix2)
		hout.write("CD1_1   =     %.8g\n" % cd1_1)
		hout.write("CD1_2   =     %.8g\n" % cd1_2)
		hout.write("CD2_1   =     %.8g\n" % cd2_1)
		hout.write("CD2_2   =     %.8g\n" % cd2_2)
		hout.write("END")
		hout.close()
		swarpcmd += '-IMAGE_SIZE "%i, %i" ' % (nax1, nax2)
    
    scmd=os.popen(swarpcmd, 'r', -1)
    scmd.readlines()
    shutil.move('coadd.fits', outfile)
    
    if rms == True:
    	swarpcmd='swarp %s.rms.fits ' % image[:-5]
    	if ractr!=None or dcctr!=None:
        	swarpcmd+='-CENTER_TYPE MANUAL -CENTER "%s, %s" ' % (ractr, dcctr)
    	if pixscale!=None:
        	swarpcmd+='-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %f ' % pixscale
    	if size!=None:
        	swarpcmd+='-IMAGE_SIZE "%i, %i" ' % (size[0], size[1])
    	if backsub==no:
        	swarpcmd+='-SUBTRACT_BACK N '
    	swarpcmd+='-COPY_KEYWORDS OBJECT,SKYSUB,SKYBKG,SEEPIX '

    	if refimage!=None:
			[nax1,nax2,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2] = get_head(refimage, ["NAXIS1", "NAXIS2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"])
			hout = open("coadd.head", "w")
			hout.write("CRVAL1  =     %.8g\n" % crval1)
			hout.write("CRVAL2  =     %.8g\n" % crval2)
			hout.write("CRPIX1  =     %.8g\n" % crpix1)
			hout.write("CRPIX2  =     %.8g\n" % crpix2)
			hout.write("CD1_1   =     %.8g\n" % cd1_1)
			hout.write("CD1_2   =     %.8g\n" % cd1_2)
			hout.write("CD2_1   =     %.8g\n" % cd2_1)
			hout.write("CD2_2   =     %.8g\n" % cd2_2)
			hout.write("END")
			hout.close()
			swarpcmd += '-IMAGE_SIZE "%i, %i" ' % (nax1, nax2)
    
    	scmd=os.popen(swarpcmd, 'r', -1)
    	scmd.readlines()
    	shutil.move('coadd.fits', '%s.rms.fits' % outfile[:-5])
    	
    
    if mask == True:
    	swarpcmd='swarp %s.mask.fits ' % image[:-5]
    	if ractr!=None or dcctr!=None:
        	swarpcmd+='-CENTER_TYPE MANUAL -CENTER "%s, %s" ' % (ractr, dcctr)
    	if pixscale!=None:
        	swarpcmd+='-PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %f ' % pixscale
    	if size!=None:
        	swarpcmd+='-IMAGE_SIZE "%i, %i" ' % (size[0], size[1])
    	if backsub==no:
        	swarpcmd+='-SUBTRACT_BACK N '
    	swarpcmd+='-COPY_KEYWORDS OBJECT,SKYSUB,SKYBKG,SEEPIX '

    	if refimage!=None:
			[nax1,nax2,crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2] = get_head(refimage, ["NAXIS1", "NAXIS2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"])
			hout = open("coadd.head", "w")
			hout.write("CRVAL1  =     %.8g\n" % crval1)
			hout.write("CRVAL2  =     %.8g\n" % crval2)
			hout.write("CRPIX1  =     %.8g\n" % crpix1)
			hout.write("CRPIX2  =     %.8g\n" % crpix2)
			hout.write("CD1_1   =     %.8g\n" % cd1_1)
			hout.write("CD1_2   =     %.8g\n" % cd1_2)
			hout.write("CD2_1   =     %.8g\n" % cd2_1)
			hout.write("CD2_2   =     %.8g\n" % cd2_2)
			hout.write("END")
			hout.close()
			swarpcmd += '-IMAGE_SIZE "%i, %i" ' % (nax1, nax2)
    
    	scmd=os.popen(swarpcmd, 'r', -1)
    	scmd.readlines()
    	
    	iraf.imexpr("a > 0 ? 1 : 0", "%s.mask.fits" % outfile[:-5],
    	            "coadd.fits")
    	os.remove("coadd.fits")
    		
    return
    	
############################################################################

def p60sdsssub(inlis, refimage, ot, distortdeg=1, scthresh1=3.0, 
               scthresh2=10.0, tu=50000, iu=50000, ig=2.3, tg=1.0, 
               stamps=None, nsx=4, nsy=4, ko=0, bgo=0, radius=10, 
               tlow=None, ilow=None, sthresh=5.0, ng=None, aperture=10.0,
               distfile=P60DISTORT):

    '''P60 Subtraction using SDSS image as reference'''

    images = iraffiles(inlis)
    images.sort()

    # Get WCS center, pixel scale, and pixel extent of reference
    [n1,n2]=get_head(refimage, ['NAXIS1','NAXIS2'])
    [[ractr,dcctr]]=impix2wcs(refimage,n1/2.0,n2/2.0)
    if not (check_head(refimage,'PIXSCALE')):
        print 'Error: Please add PIXSCALE keyword to reference image'
        return 0
    else:
        pix=get_head(refimage,'PIXSCALE')

    # Create stamps file
    if not (stamps==None):
        refstars=Starlist(stamps)
        refstars.wcs2pix(refimage)
        outf=open('ref.stamps', 'w')
        for star in refstars:
            outf.write('%.1f\t%.1f\n' % (star.xval, star.yval))
        outf.close()

    # Create OT file
    refstars=Starlist(ot)
    refstars.wcs2pix(refimage)
    outf=open('ref.coo', 'w')
    outf.write('%.2f\t%.2f\n' % (refstars[0].xval, refstars[0].yval))
    outf.close()   

    # First get good distortion correction (if succifient number of images)
#    if len(images)>5:
#        p60scampall(images, 't20', 'dist', distortdeg=3, scthresh1=3.0, 
#                    scthresh2=10.0)
#    else:
#        for image in images:
#            iraf.imcopy(image, '%s.dist.fits' % image.split('.')[0])

    for image in images:

        root=image.split('.')[0]

        # Subtract sky background if necessary
        [iskysub,iskybkg]=get_head(image,['SKYSUB','SKYBKG'])
        if (iskysub==1):
            iraf.imarith(image,'+',iskybkg,image)
            update_head(image,'SKYSUB',0)

        # First use pre-determined distortion file
        #shutil.copy(distfile, '%s.head' % root)
        #p60swarp(image, '%s.dist.fits' % root, backsub=no)

        # Run scamp
        #p60scamp('%s.dist.fits' % root, refimage=refimage, 
        p60scamp('%s.fits' % root, refimage=refimage, 
                 distortdeg=distortdeg, scthresh1=scthresh1,
                 scthresh2=scthresh2, rms=True, mask=True)

        # Run Swarp
        #p60swarp('%s.dist.fits' % root, '%s.shift.fits' % root, ractr=ractr, 
        p60swarp('%s.fits' % root, '%s.shift.fits' % root, backsub=no,
                 refimage=refimage, rms=True, mask=True)
			
        # Subtract
        if (ilow==None):
        	iraf.iterstat(image)
        	ilow = iraf.iterstat.median - 10.0 * iraf.iterstat.sigma
        if (tlow==None):
        	iraf.iterstat(refimage)
        	tlow = iraf.iterstat.median - 10.0 * iraf.iterstat.sigma
        tu = get_head(refimage, "SATURATE") * 0.90
        iu = get_head(image, "SATURATE") * 0.90
        fwhm = get_head(image, "SEEPIX")
        if fwhm > 15.0:
        	radius=15
        p60hotpants('%s.shift.fits' % root, refimage, '%s.sub.fits' % root,
                    tu=tu, iu=iu, ko=ko, bgo=bgo, nsx=nsx, 
                    nsy=nsy, radius=radius, tlow=tlow, ilow=ilow, 
                    sthresh=sthresh, ng=ng, stamps="ref.stamps", scimage=no,
                    rms=True)

        # Photometer subtracted image
        #iraf.phot('%s.sub.fits' % root, coords='ref.coo', output='%s.mag' % 
                  #root, epadu=ig, exposure='', calgorithm='none',
                  #salgorithm='median', annulus=30.0, dannulus=10.0, 
                  #weighting='constant', apertures=aperture, zmag=25.0, 
                  #interactive=no) 

    print "Exiting successfully"
    return
    
###########################################################################

def p60scamp(inlis, refimage=None, distortdeg=3, scthresh1=5.0, 
             scthresh2=10.0, match=no, cat="SDSS-R9", rms=True,
             mask=True):

    '''P60 Subtraction using scamp for alignment and HOTPANTS for
    subtraction.'''

    subimages=iraffiles(inlis)

    # Create relevant SExtractor files
    f1 = open("daofind.param", "w")
    f1.write("NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO\nFLAGS\nA_IMAGE\n")
    f1.write("B_IMAGE\nELONGATION\nFWHM_IMAGE\nCLASS_STAR\nXWIN_WORLD\n")
    f1.write("YWIN_WORLD\nERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\n")
    f1.write("ERRAWIN_WORLD\nERRBWIN_WORLD\nERRTHETAWIN_WORLD\n")
    f1.write("FLUX_AUTO\nFLUX_RADIUS\nFLUXERR_AUTO\n")
    f1.close()
    
    f2 = open("default.nnw", "w")
    nnw = """NNW

# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:       9 for profile parameters + 1 for seeing.
# outputs:      ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00
 
-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00
 1.00000e+00"""

	f2.write(nnw)
	f2.close()
    
	f3.open("default.conv", "w")
	f3.write("CONV NORM\n")
	f3.write("# 3x3 convolution mask of a gaussian PSF with FWHM = 2.0 pixels.\n")
	f3.write("0.260856 0.483068 0.260856\n")
	f3.write("0.483068 0.894573 0.483068\n")
	f3.write("0.260856 0.483068 0.260856\n")
	f3.close()
	
	f4.open("default.sex", "w")
	sx = """# Default configuration file for SExtractor V1.2
# EB 18/08/97
# (*) indicates parameters which can be omitted from this config file.

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME    test.cat        # name of the output catalog
CATALOG_TYPE    FITS_LDAC       # "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"

PARAMETERS_NAME daofind.param   # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE     CCD             # "CCD" or "PHOTO" (*)
DETECT_MINAREA  10              # minimum number of pixels above threshold
DETECT_THRESH   2.              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH 2.              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER          Y               # apply filter for detection ("Y" or "N")?
FILTER_NAME     default.conv    # name of the file containing the filter

DEBLEND_NTHRESH 32              # Number of deblending sub-thresholds
DEBLEND_MINCONT 0.005           # Minimum contrast parameter for deblending

CLEAN           Y               # Clean spurious detections? (Y or N)?
CLEAN_PARAM     1.0             # Cleaning efficiency

MASK_TYPE       CORRECT         # type of detection MASKing; can be one of
                                # NONE, BLANK, or CORRECT

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES  7               # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS 2.5, 3.5        # MAG_AUTO parameters: <Kron_fact>,<min_radius>

SATUR_LEVEL     50000.0         # level (in ADUs) at which arises saturation

MAG_ZEROPOINT   0.0             # magnitude zero-point
MAG_GAMMA       4.0             # gamma of emulsion (for photographic scans)
GAIN            2.3             # detector gain in e-/ADU.
PIXEL_SCALE     0               # size of pixel in arcsec (0=use FITS WCS info).

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM     2.0             # stellar FWHM in arcsec
STARNNW_NAME    default.nnw     # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE       64              # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE 3               # Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE  GLOBAL          # can be "GLOBAL" or "LOCAL" (*)

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE NONE            # can be one of "NONE", "BACKGROUND",
                                # "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
                                # "-OBJECTS", "SEGMENTATION", "APERTURES",
                                # or "FILTERED" (*)
CHECKIMAGE_NAME check.fits      # Filename for the check-image (*)

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK 3000            # number of objects in stack
MEMORY_PIXSTACK 300000          # number of pixels in stack
MEMORY_BUFSIZE  1024            # number of lines in buffer

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE    NORMAL          # can be "QUIET", "NORMAL" or "FULL" (*)

#------------------------------- New Stuff -----------------------------------

# Surprise!!"""
	    
	f4.write(nnw)
	f4.close()

    # Create WCS catalog from reference image (if necessary)
    if not (refimage==None):
        refroot,ext = refimage.split('.')
        if not os.path.exists('%s.cat' % refroot):
            if rms==True and os.path.exists("%s.rms.fits" % refimage[:-5]):
                os.system('sex -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE %s.rms.fits %s' % (refimage[:-5], refimage))
            else:
            	os.system('sex %s' % refimage)
            shutil.move('test.cat', '%s.cat' % refroot)

    for image in subimages:

        # Extract image root
        root=image.split('.')[0]

        # Create FITS-LDAC file from SExtractor
        if rms==True:
        	os.system('sex -WEIGHT_TYPE MAP_VAR -WEIGHT_IMAGE %s.var.fits %s' % (image[:-5], image))
        else:
	        os.system('sex %s' % image)
        shutil.move('test.cat', '%s.cat' % root)

        # Run scamp
        scampcmd="scamp %s.cat -DISTORT_DEGREES %i -SOLVE_PHOTOM N -SN_THRESHOLDS %f,%f -CHECKPLOT_DEV NULL " % (root, distortdeg, scthresh1, scthresh2)
        if refimage==None:
            scampcmd+="-ASTREF_CATALOG %s " % cat
        else:
            scampcmd+="-ASTREF_CATALOG FILE -ASTREFCENT_KEYS XWIN_WORLD,YWIN_WORLD -ASTREFCAT_NAME %s.cat -ASTREFERR_KEYS ERRAWIN_WORLD,ERRBWIN_WORLD,ERRTHETAWIN_WORLD -ASTREFMAG_KEY MAG_AUTO" % refroot
        if match:
            scampcmd+=" -MATCH Y -POSITION_MAXERR 1.0"
        else:
            scampcmd+=" -MATCH N"
        scmd=os.popen(scampcmd, 'r', -1)
        scmd.readlines()
        
        # Update header
        os.system('missfits %s' % image)
        
        # Update RMS image (if desired)
        if rms == True:
        	os.system("cp %s.head %s.rms.head" % (image[:-5], image[:-5]))
        	os.system("missfits %s.rms.fits" % image[:-5])
        	
        if mask == True:
        	os.system("cp %s.head %s.mask.head" % (image[:-5], image[:-5]))
        	os.system("missfits %s.mask.fits" % image[:-5])
        	
        os.system('rm %s.head' % root)

    print "Exiting successfully"

######################################################################

def psfphot(image, coofile, ot, wtimage="", varorder=1, clobber=globclob, 
            verbose=globver, pixtol=3.0, maxnpsf=25):

    """ perform PSF-based photometry on a single target star (SN?) at RA, Dec and  
        also on a set of comparison stars, using daophot.  simultaneously 
        perform aperture photometry on all the comparison stars (after 
        subtracting off contributions from neighbors) to enable absolute 
        photometry by comparison to aperture photometry of standard stars 
        observed in other fields """

    # Defaults / constants
    psfmult=5.0         #standard factor (multiplied by fwhm to get psfradius)
    psfmultsmall=3.0    #similar to psfmult, adjusted for nstar and substar

    # Necessary package
    iraf.imutil()
    iraf.digiphot()
    iraf.daophot()
    
    # Detect stars
    iqobjs("%s.sub.fits" % image[:-5], 1.5, 12000.0, wtimage=wtimage, skyval="0.0")
    
    root = image[:-5]
    [gain, rnoise, fwhm] = get_head(image, ["GAIN", "READN", "SEEPIX"])
    fwhm = float(fwhm); rnoise = float(rnoise)

    iraf.iterstat(image)
    
    # Saturation level
    if not check_head(image, "SATURATE"):
        saturate = 60000.0
    else:
        saturate = get_head(image, "SATURATE")
        	        
    # Update datapars and daopars
    iraf.datapars.fwhmpsf=fwhm
    iraf.datapars.sigma=iraf.iterstat.sigma
    iraf.datapars.datamin=iraf.iterstat.median-10*iraf.iterstat.sigma
    iraf.datapars.datamax=0.50*saturate
    iraf.datapars.readnoise=rnoise
    iraf.datapars.epadu=gain   
    iraf.daopars.psfrad=psfmult*fwhm
    iraf.daopars.fitrad=fwhm
    iraf.daopars.function="gauss,moffat15,moffat25,lorentz,penny1"
    iraf.daopars.varorder=varorder

    # Reference stars file
    stars = Starlist(coofile)
    stars.wcs2pix(image)
    outf = open("%s.coo.1" % image[:-5], "w")
    for star in stars:
        outf.write("%10.3f%10.3f\n" % (star.xval, star.yval))
    outf.close()

    #Aperture photometry
    iraf.daophot.phot(root,'default','default',apertures=fwhm,verify=no,
                      interac=no,verbose=verbose)

    iraf.datapars.datamax=0.50*saturate
    iraf.pstselect(root,'default','default',maxnpsf,interactive=no,
                   verify=no,verbose=verbose)

    iraf.psf(root,'default','default','default','default','default',
             interactive=no,showplots=no,verify=no,verbose=verbose)

    iraf.allstar(root,'default','default','default','default','default',
                 verify=no,verbose=verbose)
        
    # Prep for subtracted image
    iraf.iterstat("%s.sub.fits" % root)

    iraf.datapars.sigma=iraf.iterstat.sigma
    iraf.datapars.datamin=iraf.iterstat.median-10*iraf.iterstat.sigma
    iraf.datapars.datamax=saturate
                 
    # Look for source at OT location
    substars = Starlist("%s.sub.fits.stars" % image[:-5])
    otstars = Starlist(ot)
    otstars.wcs2pix("%s.sub.fits" % image[:-5])
    smatch, omatch = substars.match(otstars,tol=pixtol,useflags=no)

	# Generate coo file
    otcoo = open("%s.sub.coo.1" % image[:-5], "w")
    
    if len(smatch)==0:
    	otcoo.write("%10.3f%10.3f\n" % (otstars[0].xval, otstars[0].yval))
    else:
    	otcoo.write("%10.3f%10.3f\n" % (smatch[0].xval, smatch[0].yval))
    
    otcoo.close()
    
    iraf.daophot.phot("%s.sub.fits" % root, "%s.sub.coo.1" % image[:-5], 
                      'default', 'default', apertures=fwhm, calgorithm="none", 
                      interac=no, verify=no, verbose=verbose)
                      
    if len(smatch)==0:
    	print "No match in subtracted image: %s.sub.fits" % root
    else:
    	iraf.allstar("%s.sub.fits" % root, 'default', "%s.psf.1.fits" % root, 
                     'default', 'default', 'default', verify=no, verbose=no)

	return

###############################################################################

def ratircal(outfile, ims, refcal, filter):

	"""Calibrate RATIR imaging, returning photometry."""
	
	refstars = Starlist(refcal)
	refstars.set_mag("%sMAG" % filter.upper())
	
	outf = open(outfile, "w")
	
	for im in ims:
	
		[mjd, exp] = get_head(im, ["MJD-OBS", "EXPTIME"])
		
		# If something went wrong
		if (not os.path.exists("%s.sub.fits" % im[:-5])):
			outf.write("%10.3f%10.2f%10.3f%10.3f%10.3f\n" % (mjd, exp, 99.0, 99.0, 99.0))
		
		# Calibration depends on if OT was detected
		elif (not os.path.exists("%s.sub.fits.als.1" % im[:-5])):
		
			# Aperture photometry it is
			stars = Starlist("%s.mag.1" % im[:-5])
			refstars.wcs2pix(im)
			zp, zpu = stars.zeropt(refstars,method="mean",rejout=1)
			
			ulimf = open("%s.sub.fits.mag.1" % im[:-5])
			lines = ulimf.readlines()
			sky = float(lines[-3].strip().split()[1])
			aper = float(lines[-1].strip().split()[2])
			flux = float(lines[-1].strip().split()[3])
			exptime = float(lines[-2].strip().split()[0])
			
			ulim = -2.5 * np.log10( (np.max([flux,0.0]) + 5 * np.sqrt(aper) * sky) / exptime ) + 25.0 + zp
			outf.write("%10.3f%10.2f%10.3f%10.3f%10.3f\n" % (mjd, exp, ulim, 99.0, zpu))
			
		else:
		
			# PSF photometry
			stars = als2reg(im, 1)
			refstars.wcs2pix(im)
			zp, zpu = stars.zeropt(refstars,method="mean",rejout=1)
			grb = als2reg("%s.sub.fits.fits" % im[:-5], 1)
			outf.write("%10.3f%10.2f%10.3f%10.3f%10.3f\n" % (mjd, exp, grb[0].mag+zp, grb[0].magu, zpu))
			
	outf.close()
			
	return
	
###############################################################################

def ratir_doall(name, ra, dec, refdate, cat="SDSS-R9", varorder=1):

	"""Full pipeline for image subtractions with RATIR data."""
	
	filts = ["r", "i", "Z", "Y", "J", "H"]
	
	pixscale_dict = {'r': 0.317, 'i': 0.317, 'Z': 0.292, 'Y': 0.292, 'J': 0.292, 'H': 0.292}
	exptime_dict = {'r': 80.0, 'i': 80.0, 'Z': 67.11, 'Y': 67.11, 'J': 67.11, 'H': 67.11}
	gain_dict = {'r': 1.23, 'i': 1.23, 'Z': 2.20, 'Y': 2.20, 'J': 2.40, 'H': 2.40}
	readn_dict = {'r': 13.6, 'i': 13.6, 'Z': 14.70, 'Y': 14.70, 'J': 11.25, 'H': 11.25}
	sat_dict = {'r': 50000.0, 'i': 50000.0, 'Z': 24000.0, 'Y': 24000.0, 'J': 24000.0, 'H': 24000.0}

	dirs = glob.glob("20??????")
	dirs.remove(refdate)
	
	# Create directory structure 
	for filt in filts:
		os.mkdir(filt)
		
	# Rename "new" files
	for dir in dirs:
		if os.path.exists("%s/stack_C0_r.fits" % dir) and os.path.exists("%s/stack_C0_r.rms.fits" % dir):
			shutil.copy("%s/stack_C0_r.fits" % dir, "r/%s_r.fits" % dir)
			shutil.copy("%s/stack_C0_r.rms.fits" % dir, "r/%s_r.rms.fits" % dir)
		if os.path.exists("%s/stack_C1_i.fits" % dir) and os.path.exists("%s/stack_C1_i.rms.fits" % dir):
			shutil.copy("%s/stack_C1_i.fits" % dir, "i/%s_i.fits" % dir)
			shutil.copy("%s/stack_C1_i.rms.fits" % dir, "i/%s_i.rms.fits" % dir)
		if os.path.exists("%s/stackA_C2_ZY.fits" % dir) and os.path.exists("%s/stackA_C2_ZY.rms.fits" % dir):
			shutil.copy("%s/stackA_C2_ZY.fits" % dir, "Z/%s_Z.fits" % dir)
			shutil.copy("%s/stackA_C2_ZY.rms.fits" % dir, "Z/%s_Z.rms.fits" % dir)
		if os.path.exists("%s/stackB_C2_ZY.fits" % dir) and os.path.exists("%s/stackB_C2_ZY.rms.fits" % dir):
			shutil.copy("%s/stackB_C2_ZY.fits" % dir, "Y/%s_Y.fits" % dir)
			shutil.copy("%s/stackB_C2_ZY.rms.fits" % dir, "Y/%s_Y.rms.fits" % dir)
		if os.path.exists("%s/stackA_C3_JH.fits" % dir) and os.path.exists("%s/stackA_C3_JH.rms.fits" % dir):
			shutil.copy("%s/stackA_C3_JH.fits" % dir, "J/%s_J.fits" % dir)
			shutil.copy("%s/stackA_C3_JH.rms.fits" % dir, "J/%s_J.rms.fits" % dir)
		if os.path.exists("%s/stackB_C3_JH.fits" % dir) and os.path.exists("%s/stackB_C3_JH.rms.fits" % dir):
			shutil.copy("%s/stackB_C3_JH.fits" % dir, "H/%s_H.fits" % dir)
			shutil.copy("%s/stackB_C3_JH.rms.fits" % dir, "H/%s_H.rms.fits" % dir)

	# Rename reference files
	shutil.copy("%s/stack_C0_r.fits" % refdate, "r/ref_r.fits")
	shutil.copy("%s/stack_C0_r.rms.fits" % refdate, "r/ref_r.rms.fits")
	shutil.copy("%s/stack_C1_i.fits" % refdate, "i/ref_i.fits")
	shutil.copy("%s/stack_C1_i.rms.fits" % refdate, "i/ref_i.rms.fits")
	shutil.copy("%s/stackA_C2_ZY.fits" % refdate, "Z/ref_Z.fits")
	shutil.copy("%s/stackA_C2_ZY.rms.fits" % refdate, "Z/ref_Z.rms.fits")
	shutil.copy("%s/stackB_C2_ZY.fits" % refdate, "Y/ref_Y.fits")
	shutil.copy("%s/stackB_C2_ZY.rms.fits" % refdate, "Y/ref_Y.rms.fits")
	shutil.copy("%s/stackA_C3_JH.fits" % refdate, "J/ref_J.fits")
	shutil.copy("%s/stackA_C3_JH.rms.fits" % refdate, "J/ref_J.rms.fits")
	shutil.copy("%s/stackB_C3_JH.fits" % refdate, "H/ref_H.fits")
	shutil.copy("%s/stackB_C3_JH.rms.fits" % refdate, "H/ref_H.rms.fits")

	# Create ds9 region file for OT location
	coofile = open("Coords.reg", "w")
	coofile.write('fk5;circle(%10.5f,%10.5f,1.0") # text={%s}' % (ra, dec, name))
	coofile.close()
	
	# Get SDSS reference stars
	os.system("getsdss.pl -r 7.0 -f sdss.reg -p %10.5f %10.5f sdss.txt" % (ra, dec))
	
	# Get 2MASS stars
	os.system("getastrom.pl -d 2mass -r 7.0 -f temp.reg %10.5f %10.5f 2mass.txt" % (ra, dec))
	tmass = Starlist("temp.reg")
	for star in tmass:
		star.mags["YMAG"] = star.mags["JMAG"] + 0.5 * (star.mags["JMAG"] - star.mags["HMAG"]) + 0.08
	tmass.write("2mass.reg")
	os.remove("temp.reg")
	
	# If no SDSS coverage, copy 2mass.reg to sdss.reg (kludge until PS1)
	if os.stat("sdss.txt")[6] == 0:
		shutil.copy("ps1.reg", "sdss.reg")
	
	# Loop over filters
	for filt in filts:
		
		iraf.chdir(filt)
		
		# Update ref header
		update_head("ref_%s.fits" % filt, ["PIXSCALE", "SATURATE"], [pixscale_dict[filt], sat_dict[filt]])
		
		# Generate mask for reference
		iraf.imexpr("a == 0 ? 1 : 0", "ref_%s.mask.fits" % filt, "ref_%s.fits" % filt)
		
		# Generate variance for reference
		iraf.imcopy("ref_%s.rms.fits" % filt, "ref_%s.var.fits" % filt)
		iraf.imreplace("ref_%s.var.fits" % filt, 1.0e31, lower=INDEF, upper=0.0)
		
		# Detect objects in reference
		iqobjs("ref_%s.fits" % filt, 3.0, sat_dict[filt], wtimage="ref_%s.var.fits" % filt, skyval="0.0")
		
		# Tweak WCS in reference image
		p60scamp("ref_%s.fits" % filt, distortdeg=1, match=yes, rms=True, mask=True, cat=cat)
		os.remove("ref_%s.cat" % filt)
		
		# Get reference stars for PSF matching
		sdss = Starlist("../sdss.reg")
		ref = Starlist("ref_%s.fits.stars" % filt)
		sdss.wcs2pix("ref_%s.fits" % filt)
		a,b=sdss.match(ref, maxnum=1000)
		a.write("psf_%s.reg" % filt)
		
		# Loop over "new images"
		ims = glob.glob("????????_%s.fits" % filt)
		for im in ims:
		
			# Update headers
			nstack = int(get_head(im, "EXPTIME") / exptime_dict[filt])
			update_head(im, ["GAIN", "READN", "SATURATE"], [nstack * gain_dict[filt], readn_dict[filt] / np.sqrt(nstack), sat_dict[filt]])
			
			# Create variance files
			iraf.imcopy("%s.rms.fits" % im[:-5], "%s.var.fits" % im[:-5])
			iraf.imreplace("%s.var.fits" % im[:-5], 1.0e31, lower=INDEF, upper=0.0)
			
			# Create bad pixel mask
			iraf.imexpr("a == 0 ? 1 : 0", "%s.mask.fits" % im[:-5], im)
			
			# Update WCS for new images
			p60scamp(im, distortdeg=1, match=yes, rms=True, mask=True, cat=cat)
			os.remove("%s.cat" % im[:-5])
			
			# Detect objects and measure seeing
			iqobjs(im, 5.0, sat_dict[filt], wtimage="%s.var.fits" % im[:-5], skyval="0.0")
			
			# Match sources for PSF determination
			newstars = Starlist("%s.stars" % im)
			newstars.pix2wcs(im)
			newstars.wcs2pix("ref_%s.fits" % filt)
			c,d = newstars.match(a, maxnum=1000)
			d.write("psf_%s.%s.reg" % (filt, im[:-5]))
			
			# If there was a problem with alignment, will be no PSF stars
			if not os.stat("psf_%s.%s.reg" % (filt, im[:-5]))[6] == 0:
			
				if varorder==0:
				#if (filt=="J") or (filt=="H"):
					p60sdsssub(im, "ref_%s.fits" % filt, "../Coords.reg", stamps="psf_%s.%s.reg" % (filt, im[:-5]), nsx=2, nsy=2, ko=0, distortdeg=1)
				else:
					p60sdsssub(im, "ref_%s.fits" % filt, "../Coords.reg", stamps="psf_%s.%s.reg" % (filt, im[:-5]), nsx=5, nsy=5, ko=1, distortdeg=1)
			
		# PSF photometry
		for im in ims:
			if os.path.exists("%s.sub.fits" % im[:-5]):
				psfphot(im, "psf_%s.reg" % filt, "../Coords.reg", wtimage="ref_%s.var.fits" % filt, varorder=varorder)		

		# Photometry
		if (filt == "J") or (filt == "H") or (filt == "Y"):
			ratircal("%s.%s.dat" % (name, filt), ims, "../2mass.reg", filt)
		else:
			ratircal("%s.%s.dat" % (name, filt), ims, "../sdss.reg", filt)
			
		iraf.chdir("../")
		
	return
	