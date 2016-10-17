from __future__ import division,print_function

import numpy as np
import healpy
import matplotlib.pyplot as plt
import os

__all__ = ["SelectionFunctionTGAS","SelectionFunctionTGASRAVE"]


def lb2pix(nside, l, b, nest=True):
    """Convert l and b into pixel IDs"""
    _, __ = (np.pi/2-np.deg2rad(b), np.deg2rad(l))
    return healpy.ang2pix(nside, _, __, nest=nest)

def pix2lb(nside, pix, nest=True):
    """Convert pixel ID to l and b"""
    __, _ = healpy.pix2ang(nside,pix, nest=nest)
    return (np.rad2deg(_), np.rad2deg(np.pi/2-__))


class SelectionFunction(object):

    """Selection function base class"""

    def __init__(self, sf_str):
        sf = np.load(sf_str)
        self.map = sf['map']
        self.mag_bins = sf['mag_bins']
        self.nside = int(np.sqrt(np.shape(self.map)[0]/12.))


    def __call__(self, l, b, mag):
        """

        Evaluate the SF at (l,b) and at magnitude mag using nearest neighbour interpolation.

        Parameters
        ----------
        l: array_like
            Galactic longitude in degrees.
        b: array_like
            Galactic latitude in degrees.
        mag: array_like
            magnitude at which to evaluate SF.

        Returns
        -------

        sf: array_like
            the selection function evaluated at (l,b,mag)

        """
        #if something is off the top/bottom of the magnitude grid, then use the completeness at the inner/outer gridpoint
        mag = np.clip(mag, np.min(self.mag_bins), np.max(self.mag_bins))
        #list of the pixels inside which our points lie
        pix = lb2pix(self.nside, l, b) 
        mag_idx = np.abs(np.subtract.outer(self.mag_bins, mag)).argmin(0) #nearest neighbour points in magnitude array
        return self.map[pix, mag_idx]


    def sky_plot(self, mag):
        """

        Plot a slice through the SF at some magnitude.


        Parameters
        ----------

            mag: float
                magnitude at which to evaluate the completeness in each pixel.


        Returns
        -------

            ax: matplotlib axes object
                a mollview plot of the completeness on the sky at the chosen magnitude.

        """
        if mag > np.max(self.mag_bins) or mag < np.min(self.mag_bins):
            raise ValueError("Input magnitude outside of allowed range.")

        else:
            cmap = plt.cm.viridis
            cmap.set_under(color='white')
            idx = (np.abs(self.mag_bins - mag)).argmin()

            try:
                title = r"{0} completeness at ${1} = {2:.2f}$".format(
                    self._survey_label, self._magnitude_label, mag)

            except AttributeError:
                title = ""
            
            healpy.mollview(self.map[:,idx],
                nest=True, unit='Completeness fraction', cmap=cmap, 
                title=title, min=0., max=1.)
            return plt.gca()


    def magnitude_plot(self):
        """

        Plot the average completeness as a function of magnitude.

        Returns
        -------

            ax: matplotlib axes object
                a line plot of the average completeness as a function of magnitude.

        """
        fig, ax = plt.subplots()
        ax.plot(self.mag_bins,np.mean(self.map, axis=0))
        ax.set_xlabel(getattr(self, "_magnitude_latex_label", "Magnitude"))
        ax.set_ylabel("Average Completeness")
        ax.set_ylim((0., 1.))
        ax.set_xlim((np.min(self.mag_bins), np.max(self.mag_bins)))
        return ax


class SelectionFunctionTGAS(SelectionFunction):

    """
    Selection function for the TGAS sample as a function of 2MASS J band. 

    Parameters
    ----------

        nside (=8): int
            the resolution of the map. Must be 8 or 32.

    Returns
    -------

        sf: SelectionFunctionTGAS
            A SelectionFunctionTGAS object at the requested resolution.

    """

    _survey_label = "TGAS"
    _magnitude_label = "J"
    _magnitude_latex_label = r"$J/\mathrm{mag}$"

    def __init__(self, nside=8):
        if int(nside) not in (8, 32):
            raise ValueError("nside must be 8 or 32")

        super(SelectionFunctionTGAS, self).__init__(os.path.join(
            os.path.dirname(__file__), "maps",
            "tycho2_nside{:.0f}.dict".format(nside)))
        
        return None


class SelectionFunctionTGASRAVE(SelectionFunction):

    """
    Selection function for the TGAS-RAVE sample as a function of 2MASS J band and J-K color. 
    Parameters
    ----------
        nside (=8): int
            the resolution of the map. Must be 8 or 32.
    Returns
    -------
        sf: SelectionFunctionTGASRAVE
            A SelectionFunctionTGASRAVE object at the requested resolution.
    """

    _survey_label = "TGAS-RAVE"
    _magnitude_label = "J"
    _magnitude_latex_label = r"$J/\mathrm{mag}$"

    def __init__(self, nside=8):
        if int(nside) not in (8, 32):
            raise ValueError("nside must be 8 or 32")

        super(SelectionFunctionTGASRAVE, self).__init__(os.path.join(
            os.path.dirname(__file__), "maps",
            "tgas_rave_nside{:.0f}.dict".format(nside)))

    def __call__(self, l, b, mag, col):
        """
        Evaluate the SF at (l,b) and at magnitude mag using nearest neighbour interpolation.
        Parameters
        ----------
        l: array_like
            Galactic longitude in degrees.
        b: array_like
            Galactic latitude in degrees.
        mag: array_like
            magnitude at which to evaluate SF.
        col: array_like
            J-K color at which to evaluate SF.
        Returns
        -------
        sf: array_like
            the selection function evaluated at (l,b,mag,col)
        """
        #if something is off the top/bottom of the magnitude grid, then use the completeness at the inner/outer gridpoint
        mag = np.clip(mag, np.min(self.mag_bins), np.max(self.mag_bins))
        #list of the pixels inside which our points lie
        pix = lb2pix(self.nside, l, b) 
        mag_idx = np.abs(np.subtract.outer(self.mag_bins, mag)).argmin(0) #nearest neighbour points in magnitude array
        fracs = np.atleast_1d(self.map[pix, mag_idx])
        fracs[(np.abs(np.atleast_1d(b))<25.)&(np.atleast_1d(col)<0.5)] = 0. #set the SF to zero when the color/latitude cut enforced by RAVE isn't satisfied
        return fracs

