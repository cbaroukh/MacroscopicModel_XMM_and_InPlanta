import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from scipy.optimize import minimize



#########################
# Parameters ##to fill
#######################
filename = "data/data_XMM.xlsx"

chosen_data = ['Xylème 12_18 - Erlen 1', 'Xylème 12_18 - Erlen 2', 'Xylème 12_18 - Erlen 3']


#######################
# General parameters ###
#######################

substrates = ["Gln_Xyl", "Gluc_Xyl"]
nb_C_substrate = [5, 6]
nb_N_substrate = [2, 0]
M_substrate = [146.15, 186]

other_substrates = ['Phe', 'Val', 'Leu', 'Arg', 'Suc', 'Tyr', 'Pro', 'Lys', 'Thr', 'Iso', "Asparagine"]
nb_C_other_substrate = [9, 5, 6, 6, 12, 9, 5, 6, 4, 6, 4]
nb_N_other_substrate = [1, 1, 1, 4, 0, 1, 1, 2, 1, 1, 2]
M_other_substrate = [165.19, 117.148, 131.175, 174.204, 342.30, 181.191, 115.132, 146.190, 119.120, 131.175, 132.12]


nb_substrate = len(substrates)
nb_other_substrate = len(other_substrates)

nb_C = nb_C_substrate + nb_C_other_substrate + [4, 1] #nb carbon atoms
nb_N = nb_N_substrate + nb_N_other_substrate + [2, 0] #nb nitrogen atoms
M_N = 14 #molar mass nitrogen
M_C = 12 #molar mass carbon
M = M_substrate + M_other_substrate + [88.15, 0] #molar masses

metabolites = substrates + other_substrates + ["putr", "biom"] #metabolite names
nb_metabolites = nb_substrate + nb_other_substrate + 2 #addition de biomasse, putrescine

nb_reactions = nb_substrate #Number of reactions



#WT
params_growth_WT = np.array([(1.10898478e-02, 4.21294241e-04, 2.63210539e-01, 8.49750564e-01),
                          (4.13159097e-03, 2.42012991e-04, 1.14597371e-01, 1.20941856e-04)])

conso_coeff_other_substrate = [2.85832127e-04, 2.05239626e-04, 3.11072406e-04, 1.47043323e-04,
 7.42787597e-04, 1.39022128e-04, 9.60740956e-04, 1.86695645e-04,
 1.12618646e-04, 1.58110980e-04, 7.01713239e-04, 3.68270000e-04,
 9.17770000e-04]



#########################
# Excel import
#######################

sheet = pd.read_excel(filename, sheet_name=0, engine='openpyxl',header=None)

# find experiments
text = sheet[0].values.tolist()

exp_line_debut = [-1] + [ind for (ind, x) in enumerate(text) if isinstance(x, float) and math.isnan(x)]
exp_line_fin = exp_line_debut[1:] + [len(text)]
nb_exp = len(exp_line_debut)
exp_name = [text[x+1] for x in exp_line_debut]

# get data for each experiment
raw_data = sheet.values
a, nb_data_points = raw_data.shape

nb_SubtratesProducts = [(exp_line_fin[ind] - exp - 7)//2 for (ind, exp) in enumerate(exp_line_debut)]

temps = {}
concentrations = {}

for (ind, exp) in enumerate(exp_line_debut):
    temps[exp_name[ind]] = np.zeros((nb_metabolites, nb_data_points-1))
    concentrations[exp_name[ind]] = np.zeros((nb_metabolites, nb_data_points - 1))

    for j in range(nb_SubtratesProducts[ind]):
        Molecule = raw_data[exp + 2*j + 3][0].split()[0]
        ind_molecule = metabolites.index(Molecule)
        temps[exp_name[ind]][ind_molecule] = np.array(raw_data[exp + 2*j + 2][1:], dtype=np.float32)
        concentrations[exp_name[ind]][ind_molecule] = np.array(raw_data[exp + 2 * j + 3][1:], dtype=np.float32)

    #putrescine and biomass
    temps[exp_name[ind]][nb_metabolites - 2] = np.array(raw_data[exp + 2*j + 4][1:], dtype=np.float32)
    concentrations[exp_name[ind]][nb_metabolites - 2] = np.array(raw_data[exp + 2 * j + 5][1:], dtype=np.float32)

    temps[exp_name[ind]][nb_metabolites - 1] = np.array(raw_data[exp + 2 * j + 6][1:], dtype=np.float32)
    concentrations[exp_name[ind]][nb_metabolites - 1] = np.array(raw_data[exp + 2 * j + 8][1:], dtype=np.float32)


#########################
# Differential equation model
#######################

def MM(ksi, t, stoichMat, params):
    rates = np.zeros(nb_reactions)

    for ind in range(nb_substrate):
        if ksi[ind] <= 0:
            rates[ind] = 0.0 # to avoid negative rates and thus negative metabolite concentrations.
        else:
            rates[ind] = params[ind, 2] * ksi[ind]/(ksi[ind] + params[ind, 3]) * ksi[nb_metabolites-1]

    dksidt = stoichMat.dot(rates)

    return dksidt


#########################
# Stoichiometric matrix of the model
#######################

stoichMat = np.zeros((nb_metabolites, nb_reactions))  # stoichiometric matrix

for i in range(nb_substrate):
    stoichMat[i, i] = -params_growth_WT[i, 0]
    stoichMat[nb_metabolites - 2, i] = params_growth_WT[i, 1]
    stoichMat[nb_metabolites - 1, i] = 1

for i, other_substrate in enumerate(other_substrates):
    stoichMat[nb_substrate + i, substrates.index('Gln_Xyl')] = - conso_coeff_other_substrate[i]
    stoichMat[nb_substrate + i, substrates.index('Gluc_Xyl')] = - conso_coeff_other_substrate[i]


#########################
# Simulate data
#######################

#initialisation of error
errors = np.zeros(1)

for (ind, exp) in enumerate(exp_line_debut):

    if exp_name[ind] in chosen_data:

        #initial values
        ksi0 = np.divide(concentrations[exp_name[ind]][:, 0], nb_C)

        #time
        t = np.linspace(0, np.nanmax(temps[exp_name[ind]]) + 1, 2500)

        # solve the equation
        y = odeint(MM, ksi0, t, args=(stoichMat, params_growth_WT))  # integrate

        # remove negative values for other substrates
        for (ind_othersubstrate, other_substrate) in enumerate(other_substrates):
            for ind_time, time in enumerate(t):
                if y[ind_time, nb_substrate + ind_othersubstrate] < 0:
                    y[ind_time, nb_substrate + ind_othersubstrate] = 0


        #plot
        tot_plots = nb_SubtratesProducts[ind] + 2
        count_plot = 1

        error = np.zeros(1)


        if tot_plots > 5:
            plt.figure(figsize=(20, 20))
        else:
            plt.figure(figsize=(10, 10))

        for (ind_metabolite, metabolite) in enumerate(metabolites):

            if np.nanmean(concentrations[exp_name[ind]][ind_metabolite, :]) != 0.0:

                # compute error compared to experimental data
                y_time_metabolite = []
                for datatime in temps[exp_name[ind]][ind_metabolite, :]:
                    if np.isnan(datatime):
                        y_time_metabolite = y_time_metabolite + [float('nan')]
                    else:
                        for ind_time, time in enumerate(t):
                            if time >= datatime:
                                y_time_metabolite = y_time_metabolite + [y[ind_time, ind_metabolite]*nb_C[ind_metabolite]]
                                break


                error = np.concatenate((error, 1/np.nanmean(concentrations[exp_name[ind]][ind_metabolite, :]) * (y_time_metabolite - concentrations[exp_name[ind]][ind_metabolite, :])))

                # plot results
                plt.subplot(math.ceil(tot_plots / 2), 2, count_plot)
                count_plot = count_plot + 1
                plt.plot(t, y[:, ind_metabolite] * nb_C[ind_metabolite], 'b-')
                plt.plot(temps[exp_name[ind]][ind_metabolite, :], concentrations[exp_name[ind]][ind_metabolite, :],
                         'bo')
                plt.xlabel('time (h)')
                if metabolite == 'biom':
                    plt.ylabel(metabolite + ' (mg/L)')
                else:
                    plt.ylabel(metabolite + ' (mMC)')
                plt.title(metabolite)
                plt.grid(True)

        errors = np.concatenate((errors, error))

        plt.tight_layout()
        plt.savefig('output/png/' + exp_name[ind] + '.png')
        # plt.show()

        plt.clf()

        plt.close()


        # export to excel
        df = pd.DataFrame({'time (h)': t})

        for (ind_metabolite, metabolite) in enumerate(metabolites):
            if metabolite == 'biom':
                df['Biom (mg/L)'] = y[:, ind_metabolite]
            else:
                df[metabolite + ' (mM)'] = y[:, ind_metabolite]

        df.to_excel('output/xls/' + exp_name[ind] + '.xlsx', sheet_name= exp_name[ind], index=False)


squared_error = np.nansum(np.square(errors))

print(f'Global squared error for plot:', squared_error)
print()

