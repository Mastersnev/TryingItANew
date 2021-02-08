from sasmodels.resolution import Slit1D
from scipy.interpolate import interp1d
import numpy as np
import inspect


class SlitDimension():

    def __init__(self, empty_cell, wavelength=2.37 * 10 ** (-10), solid_source_angle=8.3 * 10 ** (-7)):
        self.I = empty_cell.I
        self.Q = empty_cell.Q
        self.omega = solid_source_angle  # sterradian
        self.labda = wavelength  # angstrom
        self.qh = self.Calc_horizontal_resolution()  # Slit horizontal dimension  1/angstrom
        self.qv = self.Calc_vertical_resolution()  # Slit vertical dimension  1/angstrom
        self.dimensions = np.array([self.qh, self.qv])

    def Calc_horizontal_resolution(self):
        func_i = interp1d(self.Q, self.I)
        q_range = np.linspace(self.Q.min(), self.Q.max(), 100000)
        i_calc = func_i(q_range)
        index_u = np.where(i_calc > self.I.max() / 2)[0][-1]
        horizontal_res = q_range[index_u] * 2  # gotta ask mr bouwman
        return horizontal_res

    def Calc_vertical_resolution(self):
        vertical_res = self.omega / ((self.labda / (2 * np.pi)) ** 2 * self.qh * 10 ** (10))
        return vertical_res * 10 ** (-10)

    def Dimensions(self):
        return self.dimensions


class PySmear():
    """
    Wrapper for pure python sasmodels resolution functions.
    """

    def __init__(self, Q, slit_dimensions, offset=None, **kwargs):
        self.Q = Q

        if 'q_range' in kwargs:
            if len(kwargs['q_range']) == 2:
                q_min = kwargs['q_range'][0]
                q_max = kwargs['q_range'][1]
            else:
                raise Exception(f"Expected two arguments received {len(kwargs['q_range'])}")
        else:
            q_min = self.Q.min()
            q_max = self.Q.max()

        self.resolution = self.slit_smear(self.Q, slit_dimensions)
        self.bin_range = self.get_bin_range(q_min, q_max)

        if offset is None:
            offset = np.searchsorted(self.resolution.q_calc, self.resolution.q[0])
        self.offset = offset

    def apply(self,I):
        """
        Apply the resolution function to the data.
        Note that this is called with iq_in matching data.x, but with
        iq_in[first_bin:last_bin] set to theory values for these bins,
        and the remainder left undefined.  The first_bin, last_bin values
        should be those returned from get_bin_range.
        The returned value is of the same length as iq_in, with the range
        first_bin:last_bin set to the resolution smeared values.
        """
        first_bin = self.bin_range[0]
        last_bin = self.bin_range[1]
        iq_in = I

        if last_bin is None: last_bin = len(iq_in)
        start, end = first_bin + self.offset, last_bin + self.offset
        q_calc = self.resolution.q_calc
        iq_calc = np.empty_like(q_calc)
        iq_calc[start:end + 1] = iq_in[first_bin:last_bin + 1]
        smeared = self.resolution.apply(iq_calc)
        return smeared

    __call__ = apply

    def get_bin_range(self, q_min, q_max):
        """
        For a given q_min, q_max, find the corresponding indices in the data.
        Returns first, last.
        Note that these are indexes into q from the data, not the q_calc
        needed by the resolution function.  Note also that these are the
        indices, not the range limits.  That is, the complete range will be
        q[first:last+1].
        """
        q = self.resolution.q
        first = np.searchsorted(q, q_min)
        last = np.searchsorted(q, q_max)
        return first, min(last, len(q) - 1)

    def slit_smear(self, Q, slit_dimensions, model=None):
        width = slit_dimensions[0]
        height = slit_dimensions[1]
        return Slit1D(Q, height, width)  # might be wrong way around
