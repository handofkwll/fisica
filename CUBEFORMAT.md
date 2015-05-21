# Datacube Format for 'Challenge' Targets
This document describes the format that should be used by those
simulating datacubes of 'challenge' objects to be observed by the 
FISICA telescope simulator.

### FITS Format
Each challenge datacube should be stored in FITS format, as a 3-d array
associated with the Primary HDU (Header/Data Unit). Additional HDUs
found in the file will be ignored. 

### Cube Size and Sampling
The spatial footprint of the model cube sets the size of the map made by the
reduction software. The footprint should be square and,
as a rule of thumb, should cover at least the central quarter of
the primary beam of the instrument at the lowest frequency observed.

To minimise artefacts the spatial sampling should be roughly one third
the Nyquist sampling for the longest baseline and highest frequency observed. 

The frequency range of the module cube need not fill the instrument band.
When the cube is read in the frequency grid of the model cube will be 
interpolated onto that measured by the instrument configuration. 
For this to work properly the frequency sampling should be higher than that 
of the instrument configuration by a factor of a few. Target fluxes
at frequencies in the instrument band not covered by the model will be set to 
zero.

The FISICA strawman instrument operates over 4 frequency bands. The following
table attaches some numbers to the rules listed above to give some idea of 
the sizes and resolutions of the model cubes required in each case.

Band | 1 | 2 | 3 | 4
-----|---|---|---|---
Wavelength Range (&mu;m)| 25-50 | 50-100 | 100-200 | 200-400
Frequency Range (cm<sup>-1</sup>)| 200-400 | 100-200 | 50-100 | 25-50
Map Size (arcsec) | 3.1 | 6.3 | 12.6 | 25.2
Sample Spacing (arcsec) | 8.6e-3 | 1.7e-2 | 3.4e-2 | 6.9e-2
Spectral Resolution | 400 | 200 | 100 | 50
Frequency Spacing (cm<sup>-1</sup>) | 0.25 | 0.25 | 0.25 | 0.25

Notes:
1. The map sizes cover half the primary beam for a 2m mirror observing
at the band minimum wavelength.
2. Spatial sample spacing is 1/3 \* Nyquist = lamba / (6 * b) where 
lambda is the band minimum wavelength and b the maximum baseline (=100m).
3. Frequency sample spacing is band minimum frequency / (2 * spectral resolution).

### An Example
A model cube simulating the pole-on view of a circumstellar disk in band 2 
was constructed by Catherine Walsh at Leiden Observatory. A summary of its 
characteristics is as follows:

No. | Name | Type | Cards | Dimensions | Format
----|------|------|-------|------------|-------
0 | PRIMARY | PrimaryHDU | 30 | (351, 351, 201) | float64   

Header information:

> SIMPLE  =                    T / file does conform to FITS standard             
> BITPIX  =                  -64 / number of bits per data pixel                  
> NAXIS   =                    3 / number of data axes                            
> NAXIS1  =                  351 / length of data axis 1                          
> NAXIS2  =                  351 / length of data axis 2                          
> NAXIS3  =                  201 / length of data axis 3                          
> EXTEND  =                    T / FITS dataset may contain extensions            
> COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
> COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H 
> CTYPE1  = 'RA---SIN'                                                            
> CTYPE2  = 'DEC--SIN'                                                            
> CTYPE3  = 'FREQ    '                                                            
> CDELT1  =        -3.462585E-06                                                  
> CDELT2  =         3.462585E-06                                                  
> CDELT3  =        9.1343302E+09                                                  
> CRPIX1  =                 176.                                                  
> CRPIX2  =                 176.                                                  
> CRPIX3  =                 100.                                                  
> CRVAL1  =                   0.                                                  
> CRVAL2  =                   0.                                                  
> CRVAL3  =        4.5640056E+12                                                  
> CUNIT1  = 'DEG     '                                                            
> CUNIT2  = 'DEG     '                                                            
> CUNIT3  = 'HZ      '                                                            
> BUNIT   = 'JY/PIXEL'                                                            
> HISTORY Image was compressed by CFITSIO using scaled integer quantization:      
> HISTORY   q = 4.000000 / quantized level scaling parameter                      
> HISTORY 'SUBTRACTIVE_DITHER_1' / Pixel Quantization Algorithm                   
> CHECKSUM= 'YkDHaiBHZiBHaiBH'   / HDU checksum updated 2013-11-12T21:06:59       
> DATASUM = '330575719'          / data unit checksum updated 2013-11-12T21:06:59 

The cube has spatial dimensions (351, 351) and spectral dimension 100. Spatial 
axes are specified in degrees from the centre, frequencies are in Hz. The 
wavelength range is 55-82&mu;m and the frequency resolution 500. The spatial
sampling is 1.25e-2 arcsec and the area covered 4.4 arcsec square. 
Examination will show these numbers to be comparable to those for band 2
in the above table.


