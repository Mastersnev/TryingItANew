import numpy as np
import scipy.interpolate as ip
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class ProcessedRawData:

    # empty cell measurement
    def __init__(self, measurement, empty_cell_measurement, true_range=True):
        """
        :param measurement: The measurement to be corrected.
        :param empty_cell_measurement: The empty cell measurement
        :param true_range: If True, it will plot the graph at minimal Q with postive Intensity.
                            Takes into account whether q_zero is bigger than the minimal postive q value.
        """
        # Making parameters easier to use
        self.ecm = empty_cell_measurement
        self.measure = measurement
        self.name = self.measure.name
        self.I_bgd = self.measure.Background()
        #self.ds = self.measure.sample_thickness #* 2.54

        temp_name = np.array(['KKB13537', 'KKB13541', 'KKB13526', 'KKB13533', 'KKB13527', 'KKB13534',
                     'KKB13524', 'KKB13531', 'KKB13529', 'KKB13535', 'KKB13530', 'KKB13536',
                     'KKB13552', 'KKB13557', 'KKB13560', 'KKB13564', 'KKB13538', 'KKB13542',
                     'KKB13559', 'KKB13563', 'KKB13558', 'KKB13562', 'KKB13550', 'KKB13556',
                     'KKB13549', 'KKB13555', 'KKB13548', 'KKB13554', 'KKB13546', 'KKB13553',
                     'KKB13545', 'KKB13551', 'KKB13539', 'KKB13543', 'KKB13540', 'KKB13544',
                     'KKB13561', 'KKB13565', ])
        temp_value = np.array([0.35, 0.35, 0.2, 0.2, 0.18, 0.18, 0.25, 0.25, 0.12, 0.12, 0.12, 0.12, 0.33, 0.33,
                      0.23, 0.23, 0.32, 0.32, 0.27, 0.27, 0.3, 0.3, 0.3, 0.3, 0.32, 0.32, 0.34, 0.34,
                      0.3, 0.3, 0.32, 0.32, 0.36, 0.36, 0.32, 0.32, 0.27, 0.27])
        #print(self.name[0:3] + self.name[5:])
        self.ds = temp_value[temp_name == self.name[0:3] + self.name[5:]]
        self.ds *= 2.54

        # Setting variables for better overview compared to using class methods.
        q_emp_0 = self.ecm.I.argmax()

        # Too lazy to find the proper q_zero, therefore I use the closest value based on max intensity or where
        # transmission is equal to zero.
        # Just a hot fix for finding q zero, To find better approx maybe interpolate or apply gauss
        if self.measure.T.argmin() >= self.measure.I.argmax():
            q_zero = self.measure.T.argmin()
        else:
            q_zero = self.measure.I.argmax()
        Iemp = self.ecm.I
        Iemp_n = ip.interp1d(self.ecm.Q, Iemp, 'linear', bounds_error=False, fill_value=Iemp[-1])

        # Making sure the q values are not out of range of the interpolation
        # Normalise intensities and discarding useless correction
        Q = self.measure.Q[self.measure.Q >= self.ecm.Q[:].min()]
        I = self.measure.I[self.measure.Q >= self.ecm.Q[:].min()]

        # Calculating the transmission data.
        # Takes average trans value of stable region, the stable region is taken to be the last 15% of the measurement
        T_emp = np.average(self.ecm.T[int(0.9 * self.ecm.T.size):])
        trans = np.average(self.measure.T[int(0.9 * len(self.measure.T)):])
        self.T_rock = self.measure.I.max() / Iemp.max()
        self.T_wide = trans / T_emp
        self.Tsas = self.T_rock / self.T_wide
        # Correcting the data
        I_corrected = I - self.T_rock * Iemp_n(Q) - (1 - self.T_rock) * self.I_bgd
        #I_corrected *= self.T_rock
        self.peak = self.ecm.I.max()

        # Finding an appropriate starting location for the graph. Negative values are discarded.
        start_array = np.where(I_corrected < 0)[0]
        if true_range:
            if start_array.size != 0:
                start = start_array[-1] + 1
                if q_zero >= start:
                    start = q_zero

            else:
                start = 0
            self.scattering_var = Q[start:]
            self.I_corrected = I_corrected[start:]
        else:
            self.scattering_var = Q
            self.I_corrected = I_corrected

    def I(self,ds,omega,peak):
        k = 1 / (ds * omega * peak * self.T_wide)
        return k * self.I_corrected

    """
    def I(self, omega=8.3e-7, absolute=True):  # 8.3e-7
        :param ds: Sample thickness in mm
        :param omega: Solid source angle in sr
        :param absolute: If True, converts to absolute units.
        :return: Returns corrected intensity of measurement

        if absolute:
            k = 1 / (self.ds * omega * self.peak * self.T_wide)
            return k * self.I_corrected

        return self.I_corrected
    """

    def Q(self):
        """
        :return: Returns Q-range used for the corrected data. Useful for plotting.
        """
        return self.scattering_var

    def T_values(self):
        """
        :return: Returns T_sas value found.
        """
        return round(self.T_rock, 3), round(self.T_wide, 3), round(self.Tsas, 3)

    def name(self):
        return self.name
