# SelectTGAS
Selection function for the TGAS sample from Gaia DR1. Requires `numpy`, `matplotlib` and `healpy`. 
The module is called `tgas_sf`.

To install using `pip`, type 

`pip install https://github.com/anguswilliams91/SelectTGAS/archive/master.zip`

The selection function is constructed by dividing the sky into equal-area bins 
using HEALPix, and then evaluating the completeness in a range of J-band magnitude bins relative to 
the 2MASS point source catalog. There are two resolutions: low (NSIDE=8) and high (NSIDE=32). 

The selection function is implemented as a class `SelectionFunctionTGAS` and the `__call__` method 
evaluates the selection function at (l,b,J) using nearest-neighbour interpolation.

Example snippet:

```python
import tgas_sf as sf

#load high-res SF for the sample
selecTGAS = sf.SelectionFunctionTGAS(nside=32)

#plot the completeness at J=8. 
selecTGAS.sky_plot(9.)

#Find the completeness at l=40 degrees , b=30 degrees , J=8
selecTGAS(40.,30.,8.)

```
![Alt text](TGAS_slice.png?raw=true)