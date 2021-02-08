import numpy as np


class RawMeasurement:

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
        self.A = array[:, 5]
        self.t = array[:, 6]
        self.name = footnotes[0].rsplit(":")[-1].replace(";", "").strip()
        self.sample = footnotes[2].rsplit(":")[-1].strip()
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

    def Name(self):
        """
        :return: Returns the name of the file.
        """
        return self.name

    def Sample(self):
        """
        :return: Name of the sample.
        """
        return self.sample

    def Sample_thickness(self):
        '''
        :return: Thickness of the sample in cm.
        '''
        return self.sample_thickness

    def Background(self):
        return self.BGD

    def Time(self):
        return self.t

    def A(self):
        """
        :return: Angle of analyser at the measurement.
        """
        return self.A


class ProcessedMeasurements:

    def __init__(self, path):
        """
        :param path: The path to the file
        """
        array = np.loadtxt(open(path, 'rt').readlines()[:-5], skiprows=4, dtype=None)
        header = np.loadtxt(open(path, 'rt').readlines()[:4], dtype=str, delimiter="\n")
        footnote = np.loadtxt(open(path, 'rt').readlines()[-5:], dtype=str, delimiter="\n")
        self.Q = array[:, 0]
        self.I = array[:, 1]
        self.name = (header[-1].rsplit(":")[1].replace(";", "")).rsplit('.')[0].strip()
        self.T_rock = float(footnote[-2].split(";")[0].split('=')[1])
        self.T_wide = float(footnote[-2].split(";")[1].split('=')[1])
        self.Tsas = float(footnote[-2].split(";")[2].split('=')[1])
        # self.sample = footnotes[2].rsplit(":")[-1]
        # self.sample_thickness = float(footnotes[4].rsplit(":")[-1])
        # self.BGD = float(footnotes[5].rsplit(":")[-1])

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
        return self.name

    def T_values(self):
        return round(self.T_rock,4) , round(self.T_wide,4), round(self.Tsas,4)