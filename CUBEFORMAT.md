# Datacube Format for 'Challenge' Targets
This document describes the format to be used for datacubes that can
be read in and then 'observed' by the FISICA telescope simulator.

The datacube should be stored in FITS format, as a 3-d array
associated with the Primary HDU (Header/Data Unit). Additional HDUs
will be ignored. 

## Cube Size and Sampling
The spatial size of the model cube sets the size of the map made by the
reduction software. The spatial footprint should be square and,
as a rule of thumb, should cover at least the central quarter of
the primary beam of the instrument at the lowest frequency observed. 
Because the primary beam size increases as frequency decreases a map 
sized to suit the lowest frequency will be good for all.

To minimise artefacts the spatial sampling should be roughly one third
the Nyquist sampling at the highest frequency observed. The Nyquist sample
spacing falls as frequency goes up so a value that works for the highest
frequency will work for all.

The frequency range of the module cube need not fill the instrument band.
The frequency sampling of the model cube will be interpolated when the
cube is read in to fit that measured by the instrument configuration. 
For this to work properly the frequency sampling should be higher than that 
of the instrument configuration by a small factor at least. Frequencies 
in the instrument band that are not specified by the model will be set to 
zero.

The FISICA strawman instrument operates over 4 frequency bands. The following
table attaches numbers to the rules listed above to give some idea of 
the sizes and resolutions of the model cubes required in each case.

Band | 1 | Band 2 | Band 3 | Band 4
Wavelength Range (um)| 25-50 | 50-100 | 100-200 | 200-400
Frequency Range (cm-1)| 200-400 | 100-200 | 50-100 | 25-50
Primary Beam Radius (arcsec) | 6.3 | 12.6 | 25.2 | 50.3
Sample Spacing (arcsec) | 8.6e-3 | 1.7e-2 | 3.4e-2 | 6.9e-2
Spectral Resolution | 400 | 200 | 100 | 50
Frequency Spacing (cm-1) | 0.25 | 0.25 | 0.25 | 0.25

Notes:
1. The primary beams are calculated assuming that each flux collector has 
a 2m primary mirror.
2. Spatial sample spacing is 1/3 * Nyquist = (lamba / 6b) where lambda is the
minimum band wavelength and b the maximu baseline (100m).
3. Frequency sample spacing is minimum frequency / (2 * spectral resolution).

## An Example
A model cube simulating the pole-on view of a circumstellar disk to be
observed in band 2 was constructed by Catherine Walsh at Leiden Observatory.  
The summary of its characteristics is as follows:

*No.    Name         Type      Cards   Dimensions   Format
*0    PRIMARY     PrimaryHDU      30   (351, 351, 201)   float64   
*SIMPLE  =                    T / file does conform to FITS standard             
*BITPIX  =                  -64 / number of bits per data pixel                  
*NAXIS   =                    3 / number of data axes                            
*NAXIS1  =                  351 / length of data axis 1                          
*NAXIS2  =                  351 / length of data axis 2                          
*NAXIS3  =                  201 / length of data axis 3                          
*EXTEND  =                    T / FITS dataset may contain extensions            
*COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
*COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H 
*CTYPE1  = 'RA---SIN'                                                            
*CTYPE2  = 'DEC--SIN'                                                            
*CTYPE3  = 'FREQ    '                                                            
*CDELT1  =        -3.462585E-06                                                  
*CDELT2  =         3.462585E-06                                                  
*CDELT3  =        9.1343302E+09                                                  
*CRPIX1  =                 176.                                                  
*CRPIX2  =                 176.                                                  
*CRPIX3  =                 100.                                                  
*CRVAL1  =                   0.                                                  
*CRVAL2  =                   0.                                                  
*CRVAL3  =        4.5640056E+12                                                  
*CUNIT1  = 'DEG     '                                                            
*CUNIT2  = 'DEG     '                                                            
*CUNIT3  = 'HZ      '                                                            
*BUNIT   = 'JY/PIXEL'                                                            
*HISTORY Image was compressed by CFITSIO using scaled integer quantization:      
*HISTORY   q = 4.000000 / quantized level scaling parameter                      
*HISTORY 'SUBTRACTIVE_DITHER_1' / Pixel Quantization Algorithm                   
*CHECKSUM= 'YkDHaiBHZiBHaiBH'   / HDU checksum updated 2013-11-12T21:06:59       
*DATASUM = '330575719'          / data unit checksum updated 2013-11-12T21:06:59 

The cube has spatial dimensions (351, 351) and spectral dimension 100. Spatial 
axes are specified in degrees from the centre, frequencies are in Hz. The 
wavelength range is 55-82um and the frequency resolution 500. The spatial
sampling is 1.25e-2 arcsec and the area covered 4.4 arcsec square.


