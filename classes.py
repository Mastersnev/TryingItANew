import numpy as np
import scipy.interpolate as ip

class Measurement:

    # noinspection PyTypeChecker
    def __init__(self, path):
        """
        :param path: The path to the file
        """
        array = np.loadtxt(open(path, 'rt').readlines()[:-8], skiprows=2, dtype=None)
        footnotes = np.loadtxt(open(path, 'rt').readlines()[-8:], dtype=str, delimiter="\n")
        self.Q = array[:, 0]
        self.I = array[:, 1]
        self.T = array[:, 3]
        self.B = array[:, 4]
        self.DT = array[:, 7].sum()
        self.filename = footnotes[0].rsplit(":")[-1].replace(";", "")
        self.sample = footnotes[2].rsplit(":")[-1]
        self.sample_thickness = float(footnotes[4].rsplit(":")[-1])
        self.BGD = float(footnotes[5].rsplit(":")[-1])

    def Q(self):
        """
        :return: Returns the array with Q-values
        """
        return self.Q

    def I(self):
        """
        :return: Returns the array with Intensities
        """
        return self.I

    def T(self):
        """
        :return: Returns the array with Transmission detector values
        """
        return self.T

    def name(self):
        """
        :return: Returns the name of the file.
        """
        return self.filename

    def sample(self):
        return self.sample

    def sample_thickness(self):
        return self.sample_thickness

    def background(self):
        return self.BGD

    def time(self):
        return self.DT

class Processed_Data:

    # empty cell measurement
    def __init__(self, measurement, empty_cell_measurement, true_range=True):
        """
        :param measurement: The measurement to be corrected.
        :param empty_cell_measurement: The empty cell measurement
        :param true_range: If True, it will plot the graph at minimal Q with postive Intensity.
                            Takes into account whether q_zero is bigger than the minimal postive q value.
        """
        # Making parameters easier to use
        ecm = empty_cell_measurement
        measure = measurement
        self.I_bgd = measure.background()
        self.ds = measure.sample_thickness

        # Setting variables for better overview compared to using class methods.
        q_emp_0 = ecm.I.argmax()


        # Too lazy to find the proper q_zero, therefore I use the closest value based on max intensity or where
        # transmission is equal to zero.
        # Just a hot fix for finding q zero, To find better approx maybe interpolate or apply gauss

        if measure.T.argmin() >= measure.I.argmax():
            q_zero = measure.T.argmin()
        else:
            q_zero = measure.I.argmax()
        Iemp = ecm.I/ecm.B/ecm.time()
        Iemp_n = ip.interp1d(ecm.Q,Iemp , 'linear')

        # Making sure the q values are not out of range of the interpolation
        # Normalise intensities and discarding useless correction

        Q1 = measure.Q[np.logical_and(measure.Q >= ecm.Q[q_emp_0:].min(), measure.Q <= ecm.Q[q_emp_0:].max())]
        Q2 = measure.Q[measure.Q > Q1.max()]
        I = measure.I / measure.B/measure.time()
        I_back = self.I_bgd / measure.B/measure.time()
        I_n1 = I[np.logical_and(measure.Q >= ecm.Q[q_emp_0:].min(), measure.Q <= ecm.Q[q_emp_0:].max())]
        I_n2 = I[measure.Q > Q1.max()]
        I_bgd_n = I_back[np.logical_and(measure.Q >= ecm.Q[q_emp_0:].min(), measure.Q <= ecm.Q[q_emp_0:].max())]
        I_bgd_n2 = I_back[measure.Q > Q1.max()]

        # Calculating the transmission data.
        # Takes average trans value of stable region, the stable region is taken to be the last 15% of the measurement

        T_emp = np.average(ecm.T[int(0.9 * ecm.T.size):]/ecm.time())
        trans = np.average(measure.T[int(0.9 * len(measure.T)):]/measure.time())
        T_rock = I.max() / Iemp.max()
        self.T_wide = trans / T_emp
        self.Tsas = T_rock / self.T_wide
        # Correcting the data
        I_corrected1 = I_n1 - T_rock * Iemp_n(Q1) - (1 - T_rock) * I_bgd_n
        I_corrected2 = I_n2 - (1 - T_rock) * I_bgd_n2
        # Finding an appropriate starting location for the graph. Negative values are discarded.
        self.peak = Iemp.max()
        start_array = np.where(I_corrected1 < 0)[0]
        if true_range:
            if start_array.size != 0:
                start = start_array[-1] + 1
                if q_zero >= start:
                    start = q_zero

            else:
                start = 0
            self.scattering_var = np.concatenate((Q1[start:],Q2))
            self.I_corrected = np.concatenate((I_corrected1[start:],I_corrected2))
        else:
            self.scattering_var = np.concatenate((Q1,Q2))
            self.I_corrected = np.concatenate((I_corrected1,I_corrected2))

    def I(self, omega=0, absolute=True):
        """
        :param ds: Sample thickness in mm
        :param omega: Solid source angle in sr
        :param absolute: If True, converts to absolute units.
        :return: Returns corrected intensity of measurement
        """
        if absolute:
            k = 1 / (self.ds * omega * self.peak * self.T_wide)

            return k*self.I_corrected

        return self.I_corrected

    def Q(self):
        """
        :return: Returns Q-range used for the corrected data. Useful for plotting.
        """
        return self.scattering_var

    def trans_sas(self):
        """
        :return: Returns T_sas value found.
        """
        return self.Tsas

class Processed_Measurements:

    def __init__(self,path):
        """
        :param path: The path to the file
        """
        array = np.loadtxt(open(path, 'rt').readlines()[:-5], skiprows=4, dtype=None)
        header = np.loadtxt(open(path, 'rt').readlines()[:4], dtype=str, delimiter="\n")
        self.Q = array[:, 0]
        self.I = array[:, 1]
        self.filename = (header[-1].rsplit(":")[1].replace(";", "")).rsplit('.')[0]
        #self.sample = footnotes[2].rsplit(":")[-1]
        #self.sample_thickness = float(footnotes[4].rsplit(":")[-1])
        #self.BGD = float(footnotes[5].rsplit(":")[-1])

    def Q(self):
        """
        :return: Returns the array with Q-values
        """
        return self.Q

    def I(self):
        """
        :return: Returns the array with Intensities
        """
        return self.I

    def name(self):
        return self.filename

    """
    def sample(self):
        return self.sample

    def sample_thickness(self):
        return self.sample_thickness

    def background(self):
        return self.BGD
"""