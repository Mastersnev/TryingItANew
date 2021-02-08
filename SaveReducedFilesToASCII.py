import numpy as np


class SaveData():

    def __init__(self, data, filelocation ):
        self.Q = data.Q
        self.I = data.I
        self.path = filelocation + '\\' +data.name + 'Reduced.txt'

    def Save_to_ascii(self):
        data_array = np.array([self.Q,self.I]).T
        np.savetxt(self.path, data_array, newline='\n',delimiter='\t')
