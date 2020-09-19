import numpy as np
import pandas as pd
import os
from scipy.stats import linregress
from sympy import diff
from sympy.abc import P


def compute_diffusivity(gas_path, mode):
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
                
                if mode == 'corrected':
                    raw_matrix = np.loadtxt(f'{root}{os.sep}{name}', usecols=[0, 1])
                    raw_matrix[:,1] **= 2
                else:
                    raw_matrix = np.loadtxt(f'{root}{os.sep}{name}', usecols=[0, 2])
                
                data.append(raw_matrix)
        
        # Se o diretório não possuir *.txt ele não entra no if
        if len(data) != 0:
            data_array = np.asarray(data)

            #data_array[:, :, 1] **= 2 # Eleva os valores da segunda coluna das matrizes ao quadrado

            filtered_array = [matrix[matrix[:,0] >= 1000] for matrix in data_array]
            
            fit_params = map(linregress, filtered_array) # Faz o ajuste linear das matrizes desprezando as primeiras linhas
            
            dc = [p.slope / 6 for p in fit_params] # Lista com as triplicatas das difusividades
            
            # Dicionário com as médias das difusividades para cada pressão
            dcm[float(key)] = np.mean(dc) * 1e-4 # A²/ps --> cm²/s

    return dict(sorted(dcm.items()))


def compute_termodynamic_factor(gas_type, press_list):
    parameters = {
        'CH4': {
            'Qm': 0.007953861011741553,
            'B': 0.027598847706793758},
        'H2': {
            'Qm': 0.01322445161143216,
            'B': 0.002778741349613243}
    }

    gas_selected = parameters[gas_type]

    qm = gas_selected['Qm']
    B = gas_selected['B']

    concentration = qm * B * P / (1 + B*P)
    concentration_diff = diff(concentration, P)

    termodynamic_factor = []

    for press in press_list:
        conc = concentration.subs(P, press)
        deriv = concentration_diff.subs(P, press)
        
        termodynamic_factor.append(float(conc / press / deriv))

    return termodynamic_factor


def compute_self_diffusivity(gases):
    """Compute gases self diffusion coefficient.

    Args:
        gases (list): path of gases folders

    Returns:
        dict: dictionary of dataframes where each one refers to a certain gas
    """
    gas_dict = None
    
    for gas in gases:
        dc_dict = compute_diffusivity(gas_path=gas, mode='self')
        gas_name = gas.split(os.sep)[-1]

        if gas_dict is None:
            gas_dict = {}

        pressures_list_in_bar = [5, 10, 25, 100]

        gas_dict[gas_name] = pd.DataFrame(data=list(dc_dict.values()), index=pressures_list_in_bar, columns=['Ds (cm^2/s)'])

    return gas_dict


def compute_transport_diffusivity(gases):
    """Compute gases transport diffusion coefficient.

    Args:
        gases (list): caminhos das pastas dos gases
        path (str): caminho da pasta da estrtura

    Returns:
        dict: Dicionário de dataframes, onde cada um é referente a um gas
    """
    gas_dict = None

    for gas in gases:
        
        dc_dict = compute_diffusivity(gas_path=gas, mode='corrected')

        gas_name = gas.split(os.sep)[-1]

        if gas_dict is None:
            gas_dict = {}

        pressures_list_in_bar = [5, 10, 25, 100]

        gas_dict[gas_name] = pd.DataFrame(data=list(dc_dict.values()), index=pressures_list_in_bar, columns=['Dc (cm^2/s)'])

        gas_dict[gas_name]['Termodynamic Factor'] = compute_termodynamic_factor(gas_type=gas_name, press_list=pressures_list_in_bar)
        
        gas_dict[gas_name]['Dt (cm^2/s)'] = gas_dict[gas_name]['Dc (cm^2/s)'] * gas_dict[gas_name]['Termodynamic Factor']

        gas_dict[gas_name] = gas_dict[gas_name].round({'Termodynamic Factor': 2})

    return gas_dict


def build_df(path, mode):

    gases = [folder.path for folder in os.scandir(path) if folder.is_dir()]

    if mode == 'corrected':
        gas_dict = compute_transport_diffusivity(gases)
    else:
        gas_dict = compute_self_diffusivity(gases)

    Diff_table = pd.concat(gas_dict, axis=1)
    Diff_table.index.name = 'P (bar)'

    return Diff_table


#t = build_df(r'C:\Users\User\Documents\MEGA\TCC\Simulações\producao\300K\1.0', mode='corrected')
#print(t)
