import numpy as np
import pandas as pd
import os
from scipy.stats import linregress
from sympy import diff
from sympy.abc import P


def compute_corrected_diffusivity(gas_path):
    """Obtém os dados dos arquivos de saída do lammps e calcula as difusividades corrigidas
    Args:
        gas_path (str): caminho da pasta do determinado gás

    Returns:
        dict[float, float]: Dicionário onde as chaves são as quantidades de moléculas e os valores são as difusividades
    """
    dcm = {}
    # Varredura recursiva em gas_path
    for root, _, files in os.walk(gas_path):
        data = []
        for name in files:
            if name.endswith('txt'):
                #print(f'{root}{os.sep}{name}')
                key = root.split(os.sep)[-2] # obtém a pressão utilizada na simulação
                # lista com as matrizes:
                # 1ª col) tempo;
                # 2ª col) Valor absoluto dos deslocamentos
                # 3ª col) Valor dos deslocamentos quadráticos
                data.append(np.loadtxt(f'{root}{os.sep}{name}', usecols=[0, 1]))
        
        # Se o diretório não possuir *.txt ele não entra no if
        if len(data) != 0:
            data_array = np.asarray(data)
            data_array[:, :, 1] **= 2 # Eleva os valores da segunda coluna das matrizes ao quadrado 
            #index = np.where(data_array[0, :, 0] == 1000)[0][0]
            fit_params = map(linregress, data_array[:, 4000:]) # Faz o ajuste linear das matrizes desprezando as primeiras linhas
            dc = [p.slope / 6 for p in fit_params] # Lista com as triplicatas das difusividades
            # Dicionário com as médias das difusividades para cada pressão
            dcm[float(key)] = np.mean(dc) * 1e-4 # A²/ps --> cm²/s

    return dict(sorted(dcm.items()))


def compute_termodynamic_factor(gas_type, press_list, temp):
    parameters = {
        'CH4': {
            'a1': 61.38, 'a2': -1.59, 'a3': 6.16e-8, 'a4': 854.66},
        'H2': {
            'a1': 1.9e-11, 'a2': 3.42, 'a3': 8.2e-12, 'a4': 2662.61}
    }

    gas_selected = parameters[gas_type]

    qm = gas_selected['a1'] * temp**gas_selected['a2']
    B = gas_selected['a3'] * np.exp(gas_selected['a4'] / temp)

    concentration = qm * B * P / (1 + B*P)
    concentration_diff = diff(concentration, P)

    termodynamic_factor = []

    for press in press_list:
        conc = concentration.subs(P, press)
        deriv = concentration_diff.subs(P, press)
        
        termodynamic_factor.append(float(conc / press / deriv))

    return termodynamic_factor

        
def compute_diffusivities(gases):
    """Compute gases transport diffusion coefficient.

    Args:
        gases (list): caminhos das pastas dos gases
        path (str): caminho da pasta da estrtura

    Returns:
        dict: Dicionário de dataframes, onde cada um é referente a um gas
    """
    gas_dict = None

    for gas in gases:
        
        dc_dict = compute_corrected_diffusivity(gas_path=gas)

        gas_name = gas.split(os.sep)[-1]

        if gas_dict is None:
            gas_dict = {}

        pressures_list_in_bar = [5, 10, 25]
        convert_to_Pa = list(map(lambda x: x*1e5, pressures_list_in_bar))

        gas_dict[gas_name] = pd.DataFrame(data=list(dc_dict.values()), index=pressures_list_in_bar, columns=['Dc (cm^2/s)'])

        gas_dict[gas_name]['Termodynamic Factor'] = compute_termodynamic_factor(gas_type=gas_name, press_list=convert_to_Pa, temp=300)
        
        gas_dict[gas_name]['Dt (cm^2/s)'] = gas_dict[gas_name]['Dc (cm^2/s)'] * gas_dict[gas_name]['Termodynamic Factor']

        gas_dict[gas_name] = gas_dict[gas_name].round({'Termodynamic Factor': 2, 'Dt (cm^2/s)': 5})

    return gas_dict


def build_df(path):

    gases = [folder.path for folder in os.scandir(path) if folder.is_dir()]

    gas_dict = compute_diffusivities(gases)

    Diff_table = pd.concat(gas_dict, axis=1)
    Diff_table.index.name = 'P (bar)'

    return Diff_table


t = build_df(r'C:\Users\user\Documents\MEGA\TCC\Simulações\producao\300K\1.0')
print(t)
