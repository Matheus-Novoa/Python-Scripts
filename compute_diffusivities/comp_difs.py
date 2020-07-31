import numpy as np
import pandas as pd
import os
from scipy.stats import linregress
from sympy import diff
from sympy.abc import P
import json

def ComputeCorrectedDiffusivity(gas_path):
    dcm = {}
    for root, _, files in os.walk(gas_path):
        data = []
        for name in files:
            if name.endswith('txt'):
                print(f'{root}{os.sep}{name}')
                key = root.split(os.sep)[-2]
                data.append(np.loadtxt(f'{root}{os.sep}{name}', usecols=[0, 1]))

        if len(data) != 0:
            data_array = np.asarray(data)
            data_array[:, :, 1] **= 2
            #index = np.where(data_array[0, :, 0] == 1000)[0][0]
            fit_params = map(linregress, data_array[:, 4000:])
            dc = [p.slope / 6 for p in fit_params]
            dcm[float(key)] = np.mean(dc) * 1e-4 # A²/ps --> cm²/s

    return dict(sorted(dcm.items()))

        
def compute_diff(gases, path):
    """
    Compute gases transport diffusion coefficient.
    """

    structure = path.split(os.sep)[-1]
    Iso_data = pd.read_excel('Dados_Isotermas.xlsm', sheet_name=structure, index_col=0)
    Iso_data.reset_index(inplace=True)

    with open('parametros_isotermas.json') as f:
        IsoParams = json.load(f)

    StructureParams = IsoParams[structure]

    gas_dict = None

    for gas in gases:
        
        dc_dict = ComputeCorrectedDiffusivity(gas)

        gas_name = gas.split(os.sep)[-1]

        A = StructureParams[gas_name]['A']
        B = StructureParams[gas_name]['B']

        Iso = A * B * P / (1 + B * P)
        Iso_diff = diff(Iso, P)

        if gas_dict is None:
            gas_dict = {}

        gas_dict[gas_name] = pd.DataFrame(list(dc_dict.values()), index=dc_dict.keys(), columns=['D0 (cm^2/s)'])
        gas_dict[gas_name]['Termodyamic Factor'] = [round(float(C / Press / Iso_diff.subs(P, Press)), 2) for Press, C in
                                                    zip(Iso_data['P [bar]'], Iso_data[gas_name])]
        gas_dict[gas_name]['Dt (cm^2/s)'] = round(
            gas_dict[gas_name]['D0 (cm^2/s)'] * gas_dict[gas_name]['Termodyamic Factor'], 5)

    return gas_dict


def build_df(path):

    gases = [folder.path for folder in os.scandir(path) if folder.is_dir()]

    gas_dict = compute_diff(gases, path)

    Dc_table = pd.concat(gas_dict, axis=1)
    Dc_table.index.name = 'P (bar)'

    return Dc_table

#df = build_df(r'C:\Users\User\Dropbox_Matheus\Dropbox\Difusividade\c48a')
#print(df)
