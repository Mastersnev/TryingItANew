from ReadingData import *
from RawDataAnalysis import *
from GraphPlotter import *
from SlitSmearing import *
from SaveReducedFilesToASCII import *
from GuinierPorod import *
import numpy as np
import os
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from math import sqrt, inf
import matplotlib.pyplot as plt
"""
Structure of raw measurement array:
q, filename, y_error, data_trans, beam_monitor, m2om, det_time, actual_time, diff_time
"""


def main(**kwargs):
    """
    Handle kwargs
    """
    if 'plot' in kwargs:
        if kwargs['plot']:
            plot = True
        else:
            plot = False
    else:
        plot = False

    if 'SAS_save' in kwargs:
        if kwargs['SAS_save']:
            sas_save = True
        else:
            sas_save = False
    else:
        sas_save = False

    """
    Variables for setup
    """
    file_loc = "C:\\Users\\svene\\OneDrive\\Minor\\7.Project\\Measurements"
    save_loc = "C:\\Users\\svene\\OneDrive\\Minor\\7.Project\\Measurements\\3 - Reduced Igor Pro Files"
    sub_folder_raw = "2 - As Measured Data"
    sub_folder_pro = "4 - Absolute Intensity Scaled Data CUT"

    """
    Reading and treating the data files
    """
    data_raw = np.array(os.listdir(file_loc + '\\' + sub_folder_raw))  # finds the names of files in folder
    data_pro = np.array(os.listdir(file_loc + '\\' + sub_folder_pro))
    as_measured = np.empty(data_raw.shape, dtype=np.object)  # creates empty array for object storage
    processed = np.empty(data_pro.shape, dtype=np.object)

    for i in range(len(data_raw)):  # Storing the objects in arrays
        as_measured[i] = RawMeasurement(file_loc + '\\' + sub_folder_raw + '\\' + data_raw[i])
    for i in range(len(data_pro)):
        processed[i] = ProcessedMeasurements(file_loc + '\\' + sub_folder_pro + '\\' + data_pro[i])

    """
    Save Reduced data for use in SASview
    """
    if sas_save:
        for data in processed:
            SaveData(data, save_loc).Save_to_ascii()

    """
    Find Empty Cell Measurement
    """
    eci = -1  # empty cell index
    for i in range(len(as_measured)):
        if as_measured[i].sample == "Empty Cell":
            eci = i

    if eci == -1:
        raise Exception('Empty Cell could not be found')

    """
    Convert opened data into corrected data
    """

    results = np.empty(as_measured.size - 1, dtype=np.object)  # creates empty array for object storage
    j = 0
    empty_cell = as_measured[eci]  # defines empty cell measurement since this one isn't corrected
    for measurement in as_measured:
        if measurement != empty_cell:
            results[j] = ProcessedRawData(measurement, empty_cell, true_range=True)
            j += 1


    def I(Q,ds,omega,peak):

        return results[0].I(ds,omega,peak)

    cv = curve_fit(I,processed[0].Q,processed[0].I)[0]
    print(cv)
    Plotlog((processed[0].Q,processed[0].Q),(processed[0].I, I(processed[0].Q,cv[0],cv[1],cv[2])))
    """
    Plot the results in log scale graphs
    """
    if plot:
        for i in range(len(results)):

            if processed[i].name != results[i].name:
                print('The following graphs do not depict the same measurement')
                print(f'Processed data file :{results[i].name} \nIgor Pro file: {processed[i].name}')


            print(results[i].T_values(), processed[i].T_values())
            Plotlog((results[i].Q(), processed[i].Q), (results[i].I(), processed[i].I), title=results[i].name,
                    x_limit=[3.5 * 10 ** (-5), 10 ** (-2)], x_label='Q[$Å^{-1}$]',
                    y_label='Intensity(Q)[$cm^{-1}sr^{-1}$]',
                    labels=['Verification', 'Igor Pro Routine'])


if __name__ == "__main__":
    main(plot=True)

"""

def RD(y,delta):

    if abs(y) <= 1:
        rd = 1
    else:
        rd = (abs(y) - sqrt(y ** 2 - 1))**2


    if abs(y-delta) <= 1:
        rd1 = 1 - delta
    else:
        rd1 = (abs(y - delta) - sqrt((y - delta) ** 2 - 1)) ** 2

    return (rd*rd1)**5



q = measurements[1].Q #[np.abs(measurements[1].Q) < 2.9 * 10 ** (-5)]
i = measurements[1].I #[np.abs(measurements[1].Q) < 2.9 * 10 ** (-5)]
aoi = q * 4.74 * 10 ** (-10) / 4 / np.pi
theta = -25
q_zero = 179.55774

r = int.quad(RD,-inf,inf,args=(2.09,))[0]

print(r)
plt.plot(q,i)
"""

    #def calc_h(Q, length, radius,gyration, porod, ps, pf):
    #    a = 0.01
    #    sd = np.array([3.06251625e-05, 0.0586])
    #    I = SANSGP(Q, a, length, radius, gyration, porod, ps, pf).I
    #    I_s = PySmear(Q, sd).apply(I)
    #    return I_s
#
    #def calc_h2(Q, length, radius, porod, G, G1):
    #    # a = np.pi / 2 + 0.174533
    #    sd = np.array([3.06251625e-05, 1.90485856e-01])
    #    I = USANSGP(length=length, radius=radius, porod=porod, G=G, G_1=G1).Powerlaw(Q)
    #    I_s = PySmear(Q, sd).apply(I)
    #    return I_s
#
    #a, b = (6, 6)
    ##param, I_s1 = curve_fit(calc_h, processed[a].Q[processed[a].Q>1e-3], processed[a].I[processed[a].Q>1e-3], maxfev=5000)
    #param, I_s1 = curve_fit(calc_h, processed[a].Q, processed[a].I,maxfev=5000)
    #print(param)
    #Plotlog((processed[a].Q, processed[a].Q),
    #        (processed[a].I, calc_h(processed[a].Q, param[0], param[1], param[2], param[3], param[4],param[5])),
    #        labels=['original','smeared'])
