import numpy as np


class USANSGP():

    def __init__(self, Type="Cylinder", **kwargs):

        """

        :param Type:
        :param kwargs:
        """

        """
        Find geometry specific parameters
        """
        self.cut_off = 10 ** (-5)
        if Type == "Cylinder":
            if 'radius' in kwargs:
                self.r = kwargs['radius']
            else:
                raise Exception("Missing geometry: radius")

            if 'length' in kwargs:
                self.l = kwargs["length"]
            else:
                raise Exception("Missing geometry: length")

            if 'G' in kwargs:
                self.G = kwargs["G"]
            else:
                raise Exception("Missing factor: G")

            if 'G_1' in kwargs:
                self.G_1 = kwargs["G_1"]
            else:
                raise Exception("Missing factor: G1")
            if 'porod' in kwargs:
                self.d = kwargs["porod"]
            else:
                raise Exception("Missing exponent: Porod")

            self.s, self.s1, self.s2 = (1, 1, 0)
            self.rg1, self.rg2 = self.Radius_gyration_cylinder()
            self.rg = self.rg1  # for Q1
            if self.d<=self.s:
                self.d = self.s + 0.1

            # self.d = 1   Not sure yet
        '''
        Find transition points
        '''
        self.Q_1 = self.Calc_Q1()
        self.Q_2 = self.Calc_Q2()

    def Radius_gyration_cylinder(self):
        rg1 = self.r / (np.sqrt(2))
        rg2 = np.sqrt(self.l ** 2 / 12 + self.r ** 2 / 2)
        return rg1, rg2

    def Calc_Q1(self):

        Q_1 = 1 / self.rg * np.sqrt((self.d - self.s) * (3 - self.s) / 2)
        return Q_1

    def Calc_Q2(self):
        Q_2 = np.sqrt((self.s1 - self.s2) / (2 / (3 - self.s2) * self.rg2 ** 2 - 2 / (3 - self.s1) * self.rg1 ** 2))

        return Q_2

    def Calc_intensity(self, Q):
        I_1 = (self.G_1 / (Q ** (self.s1)) * np.exp(-(Q ** 2 * self.rg1 ** 2) / (3 - self.s1)))[
            np.logical_and(Q >= self.Q_2, Q <= self.Q_1)]

        G_2 = self.G_1 * np.exp(
            -(self.Q_2 ** 2) * (self.rg1 ** 2 / (3 - self.s1) - self.rg2 ** 2 / (3 - self.s2))) * self.Q_2 ** (
                      self.s2 - self.s1)

        I_2 = (G_2 / (Q ** (self.s2)) * np.exp(-(Q ** 2 * self.rg2 ** 2) / (3 - self.s2)))[Q <= self.Q_2]

        D = self.G * np.exp(-(self.Q_1 ** 2 * self.rg ** 2) / (3 - self.s)) * self.Q_1 ** (self.d - self.s)

        I_3 = (D / (Q ** self.d))[Q >= self.Q_1]

        I_total = np.concatenate((I_2, I_1, I_3))
        Q_range = np.concatenate((Q[Q <= self.Q_2], Q[np.logical_and(Q >= self.Q_2, Q <= self.Q_1)], Q[Q >= self.Q_1]))
        return I_total, Q_range

    def Powerlaw(self, Q):

        D = self.G * np.exp(-(self.Q_1 ** 2 * self.rg ** 2) / (3 - self.s)) * self.Q_1 ** (self.d - self.s)
        I_1 = D / (Q ** self.d)
        I = I_1

        return I


class SANSGP():

    def __init__(self, Q, alpha, length, radius,gyration, porod, ps, pf):
        self.a = alpha
        self.l = length
        self.r = radius
        self.d = porod
        self.q = Q
        self.rg = gyration
        bkgd = 0.18

        vf = np.pi * 2 * self.r ** 2 * self.l / 4
        vs = 4 * np.pi * self.rg ** 3 / 3
        fl = self.lengthfactor()
        fr = self.radiusfactor()
        fs = self.spherefactor()
        ifibre = pf * vf * fl * fr
        isphere = ps * vs * fs
        self.I = ifibre + isphere + bkgd

    def lengthfactor(self):
        # Define Variables
        fl = np.zeros(np.shape(self.q))
        # Do Calculations
        transpoint = np.sqrt(12) / (self.l * np.cos(self.a))
        indexl = np.where(self.q <= transpoint)[0]
        indexu = np.where(self.q > transpoint)[0]
        fl[indexl] = np.exp(-(self.q[indexl] * np.cos(self.a) * self.l) ** 2 / 12)
        fl[indexu] = (np.sqrt(12 / np.e) / (self.q[indexu] * np.cos(self.a) * self.l)) ** 2
        return fl

    def radiusfactor(self):
        fr = np.zeros(np.shape(self.q))
        # Do Calculations
        transpoint = np.sqrt(18) / (2 * self.r * np.sin(self.a))
        indexl = np.where(self.q <= transpoint)[0]
        indexu = np.where(self.q > transpoint)[0]
        fr[indexl] = np.exp(-(self.q[indexl] * np.sin(self.a) * 2 * self.r) ** 2 / 12)
        fr[indexu] = (np.sqrt(18 / np.e) / (self.q[indexu] * np.sin(self.a) * 2 * self.r)) ** 3
        return fr

    def spherefactor(self):
        # Define Variables
        fs = np.zeros(np.shape(self.q))
        # Do Calculations
        q1 = (1 / self.rg) * np.sqrt((3 * self.d) / 2)
        transpoint = q1
        indexl = np.where(self.q <= transpoint)[0]
        indexu = np.where(self.q > transpoint)[0]
        fs[indexl] = np.exp(-(self.q[indexl] * self.rg) ** 2 / 3)
        fs[indexu] = np.exp(-(q1 * self.rg) ** 2 / 3) * q1 ** self.d / (self.q[indexu] ** self.d)
        return fs

    def I(self):
        return self.I


"""
    def Powerlaw(self, Q, G, d, s):
        
        D = G * np.exp(-(Q1 ** 2 * self.rg ** 2) / (3 - s)) * Q1 ** (d - s)
        I_1 = D / (Q ** d)
        I = I_1
        Q_range = Q
        return I, Q_range, Q1
"""
