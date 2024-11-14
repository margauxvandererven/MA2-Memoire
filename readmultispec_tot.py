from astropy.io import fits as pyfits
from astropy import units as u
import numpy as np

def readmultispec_tot(fitsfile, reform=True, quiet=False):
    """Read IRAF echelle spectrum in multispec format from a FITS file
    
    Can read most multispec formats including linear, log, cubic spline,
    Chebyshev or Legendre dispersion spectra
    
    If reform is true, a single spectrum dimensioned 4,1,NWAVE is returned
    as 4,NWAVE (this is the default.)  If reform is false, it is returned as
    a 3-D array.
    """

    fh = pyfits.open(fitsfile)
    try:
        header = fh[0].header
        flux_ = fh[0].data
        # mean = sum(flux_)/len(flux_)
        flux = flux_/max(flux_)
    finally:
        fh.close()
    temp = flux.shape
    nwave = temp[-1]
    if len(temp) == 1:
        nspec = 1
    else:
        nspec = temp[-2]

    # get wavelength parameters from multispec keywords
    try:
        wat2 = header['wat2****']
        count = len(wat2)
    except KeyError:
        raise ValueError('Cannot decipher header, need either WAT2_ or CRVAL keywords')

    # concatenate them all together into one big string
    watstr = []
    for i in range(len(wat2)):
        # hack to fix the fact that older pyfits versions (< 3.1)
        # strip trailing blanks from string values in an apparently
        # irrecoverable way
        # v = wat2[i].value
        v = wat2[i]
        v = v + (" " * (68 - len(v)))  # restore trailing blanks
        watstr.append(v)
    watstr = ''.join(watstr)
    # find all the spec#="..." strings
    specstr = [''] * nspec

    for i in range(nspec):
        sname = 'spec' + str(i + 1)
        p1 = watstr.find(sname)
        p2_ = watstr.find('INDEF', p1)
        p2 = watstr.find('INDEF', p2_+1)
        p3 = watstr.find('"', p2 + 1)
        if p1 < 0 or p2 < 0 or p3 < 0:
            print(p1, p2, p3)
            raise ValueError('Cannot find ' + sname + ' in WAT2_* keyword')
        specstr[i] = watstr[p2+5:p3]
    
    w1 = np.asarray(specstr[i].split(), dtype=float)
    wavelen = w1[4:]
    return {'flux': flux, 'wavelen': wavelen, 'header': header}