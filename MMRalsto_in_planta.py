import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

from scipy.integrate import odeint
from scipy.interpolate import interp1d

#########################
# Parameters ##to fill
#######################
filename = "data/data_gerlin_et_al_2021.xlsx"
datasheet_values = "Cinetique_infected"
datasheet_values_std = "Cinetique_infected_std"
datasheet_transpi = "cinetique_transpi"
datasheet_fw = "cinetique_fw"


t_end = 7 #day
t_step = 2500 #nb steps

mortality = 8.52472742e-03 #0
dilution_putr = 1.73017153e-01 #0



#########################
# Excel import
#######################

xls = pd.ExcelFile(filename)
data_exp = xls.parse(datasheet_values, header=0)
data_exp_std = xls.parse(datasheet_values_std, header=0)
data_transpi = xls.parse(datasheet_transpi, header=0)
data_fw = xls.parse(datasheet_fw, header=0)


#######################
# General parameters ###
#######################

substrates = ["Gln_Xyl", "Gluc_Xyl"]
nb_C_substrate = [5, 6]
concentrations_substrates = [3.29, 0.010]

params_growth = np.array([(1.10898478e-02, 4.21294241e-04, 2.63210539e-01, 8.49750564e-01),
                          (4.13159097e-03, 2.42012991e-04, 1.14597371e-01, 1.20941856e-04)])


other_substrates = ['Phe', 'Val', 'Leu', 'Arg', 'Suc', 'Tyr', 'Pro', 'Lys', 'Thr', 'Iso', "Asparagine", "Aspartate"]

nb_C_other_substrate = [9, 5, 6, 6, 12, 9, 5, 6, 4, 6, 4, 4]

conso_coeff_other_substrate = [2.85832127e-04, 2.05239626e-04, 3.11072406e-04, 1.47043323e-04,
 7.42787597e-04, 1.39022128e-04, 9.60740956e-04, 1.86695645e-04,
 1.12618646e-04, 1.58110980e-04, 7.01713239e-04, 3.68270000e-04]


densite_plante = 0.8 #g/mL
percentage_xyl = 1.063/100  #percentage surface/surface

convCFU_mLXyl__CFU_gFW = densite_plante/percentage_xyl #CFU/mLXyl
convDO__CFU_mL = 9.17E-10 #DO
convmg_L__DO = 414 #mg/L

convDO__CFU_gFW = convCFU_mLXyl__CFU_gFW * convDO__CFU_mL
#print("convDO__CFU_gFW", convDO__CFU_gFW)
convmg_L__CFU_gFW = convCFU_mLXyl__CFU_gFW*convDO__CFU_mL*convmg_L__DO
#print("convmg_L__CFU_gFW", convmg_L__CFU_gFW)
convmg_L_CFU_mL = convDO__CFU_mL * convmg_L__DO
#print("convmg_L_CFU_mL", convmg_L_CFU_mL)



#######################
# Print kinetic parameters
#######################

nb_substrate = len(substrates)
nb_other_substrate = len(other_substrates)

for i in range(nb_substrate):
    print('Substrate: ', substrates[i])
    print('Mu_max (h-1): ', params_growth[i, 2])
    print('Mu_max (d-1): ', params_growth[i, 2]*24)
    print('Ks (mM): ', params_growth[i, 3])
    print('YS/B (mmol/gB): ', params_growth[i, 1] * 1000)
    print('Yputr/B (mmol/gB): ', params_growth[i, 1] * 1000)

    print('Flux substrate (mmol.gB-1.h-1): ', params_growth[i, 0] * params_growth[i, 2] * concentrations_substrates[i]/(params_growth[i, 3] + concentrations_substrates[i]) * 1000)
    print('Flux substrate (mmol.gB-1.d-1): ', params_growth[i, 0] * params_growth[i, 2] * concentrations_substrates[i]/(params_growth[i, 3] + concentrations_substrates[i]) * 1000 * 24)
    print('Flux putrescine (mmol.gB-1.h-1): ', params_growth[i, 1] * params_growth[i, 2] * concentrations_substrates[i]/(params_growth[i, 3] + concentrations_substrates[i]) * 1000)
    print('Flux putrescine (mmol.gB-1.d-1): ', params_growth[i, 1] * params_growth[i, 2] * concentrations_substrates[i]/(params_growth[i, 3] + concentrations_substrates[i]) * 1000 * 24)
    print('Flux_biomass (h-1): ', params_growth[i, 2] * concentrations_substrates[i] / (params_growth[i, 3] + concentrations_substrates[i]))
    print('Flux_biomass (d-1): ', params_growth[i, 2] * concentrations_substrates[i] / (params_growth[i, 3] + concentrations_substrates[i]) * 24)
    print()



ind_Gln_Xyl = substrates.index('Gln_Xyl')
ind_Gluc_Xyl = substrates.index('Gluc_Xyl')

for i in range(nb_other_substrate):
    print('Complementary substrate: ', other_substrates[i])
    print('Coeff conso: ', conso_coeff_other_substrate[i])
    print('Flux conso (mM.gB-1.h-1): ', conso_coeff_other_substrate[i] *
          (params_growth[ind_Gln_Xyl, 2] * concentrations_substrates[ind_Gln_Xyl] / (params_growth[ind_Gln_Xyl, 3] + concentrations_substrates[ind_Gln_Xyl]) + params_growth[ind_Gluc_Xyl, 2] * concentrations_substrates[ind_Gluc_Xyl] / (params_growth[ind_Gluc_Xyl, 3] + concentrations_substrates[ind_Gluc_Xyl]) * 1000))
    print()



#########################
# Metabolic network
#######################

nb_C = nb_C_substrate + nb_C_other_substrate + [4, 4, 1] #nb carbon atoms

metabolites = substrates + other_substrates + ["3HB", "putr", "biom"] #metabolite names
nb_metabolites = len(metabolites) #addition of biomasss and putrescine

nb_reactions = nb_substrate #Number of chemical reactions

stoichMat = np.zeros((nb_metabolites, nb_reactions)) #stoichiometric matrix

for i in range(nb_substrate):
    stoichMat[i, i] = -params_growth[i, 0]
    stoichMat[nb_metabolites-2, i] = params_growth[i, 1]
    stoichMat[nb_metabolites-1, i] = 1

for i, other_substrate in enumerate(other_substrates):
    stoichMat[nb_substrate + i, substrates.index('Gln_Xyl')] = - conso_coeff_other_substrate[i]
    stoichMat[nb_substrate + i, substrates.index('Gluc_Xyl')] = - conso_coeff_other_substrate[i]

print(stoichMat)



#########################
# Differential equation model
#######################

#interpolation of transpiration rate and fresh weight necessary for simulation of odes
interp_func_transpi = interp1d(data_transpi['dpi'], data_transpi['transpi_tot (mL)'])
inter_func_fwplante = interp1d(data_fw['Time (dpi)'], data_fw['fwplante (g)'])


#initial conditions

biom_init = 10**(5.66)# CFU/gFW
t_start = 3 #day

#estimation of t_start
#t_start = 0 #day
#biom_init = 1/inter_func_fwplante(1)# CFU/gFW

#t_start = 0 #day
#biom_init = 10/inter_func_fwplante(1)# CFU/gFW

#t_start = 0 #day
#biom_init = 21/inter_func_fwplante(1)# CFU/gFW

#t_start = 0.31 #day
#biom_init = 100/inter_func_fwplante(1)# CFU/gFW

#t_start = 0.669 #day
#biom_init = 500/inter_func_fwplante(1)# CFU/gFW

#t_start = 0.83 #day
#biom_init = 1000/inter_func_fwplante(1)# CFU/gFW


#mortality at 0
#t_start = 0.95 #day
#biom_init = 1/inter_func_fwplante(1)# CFU/gFW

#t_start = 1.245 #day
#biom_init = 10/inter_func_fwplante(t_start)# CFU/gFW

#t_start = 1.337 #day
#biom_init = 21/inter_func_fwplante(t_start)# CFU/gFW

#t_start = 1.535 #day
#biom_init = 100/inter_func_fwplante(t_start)# CFU/gFW

#t_start = 1.737 #day
#biom_init = 500/inter_func_fwplante(t_start)# CFU/gFW

#t_start = 1.827 #day
#biom_init = 1000/inter_func_fwplante(t_start)# CFU/gFW




print()
print('Initial biomass (cfU/gFW): ', biom_init)
print('Initial biomass (mg/L): ', biom_init * convmg_L__CFU_gFW)


#storage variables for interpolation data
ts = []
vxyls = []
transpis = []
fws = []
Ds = []


#ODEs
def MM(ksi, t):

    rates = np.zeros(nb_reactions)
    ind_biom = metabolites.index('biom')
    ind_putr = metabolites.index('putr')
    ind_3HB = metabolites.index('3HB')

    #compute rates of chemical reactions
    for ind in range(nb_substrate):
        rates[ind] = params_growth[ind, 2] * ksi[ind]/(ksi[ind] + params_growth[ind, 3]) * ksi[ind_biom]

    #compute dilution rate
    try:
        transpi_plante = interp_func_transpi(t/24) #mL/day
    except:
        if t /24 < 1:
            transpi_plante = interp_func_transpi(1 * 24/24)
        elif t / 24 > 7:
            transpi_plante = interp_func_transpi(7 * 24 / 24)  #mL/day
    try:
        fw_plante = inter_func_fwplante(t/24) #g FW
    except:
        if t / 24 < 1:
            fw_plante = inter_func_fwplante(1 * 24 / 24) #g FW
        elif t / 24 > 7:
            fw_plante = inter_func_fwplante(7 *24 / 24) #g FW



    V_plante = fw_plante / densite_plante  # mL
    V_xyl = percentage_xyl * V_plante #mL
    dilution_rate = (transpi_plante/24)/(V_xyl) #h-1

    ts.append(t/24)
    transpis.append(transpi_plante)
    fws.append(fw_plante)
    vxyls.append(V_xyl)
    Ds.append(dilution_rate)

    #odes for all metabolites
    dksidt = stoichMat.dot(rates) - dilution_rate * ksi + dilution_rate * ksiin

    #modification for biomass and putrescine (mortality and uptake of putrescine by the plant)
    dksidt[ind_biom] = stoichMat[ind_biom, :].dot(rates) - mortality * dilution_rate * ksi[ind_biom]
    dksidt[ind_putr] = stoichMat[ind_putr, :].dot(rates) - dilution_putr * dilution_rate * ksi[ind_putr]

    #modification for 3HB
    dksidt[ind_3HB] = 23/100 * mortality * dilution_rate * ksi[ind_biom]

    return dksidt


#########################
# Simulate data
#######################
print()

#initial values
ksi0 = np.zeros((nb_metabolites))

#incoming metabolites
ksiin = np.zeros((nb_metabolites))

for metab in metabolites:

    if metab in data_exp.keys():
        ind_molecule = metabolites.index(metab)
        ksi0[ind_molecule] = data_exp[metab][0]/nb_C[ind_molecule]
        ksiin[ind_molecule] = data_exp[metab][0]/nb_C[ind_molecule]

        if metab != 'biom' and metab != 'putr':
            print('Flux entrant ', metab, ' :', ksiin[ind_molecule])


ksiin[metabolites.index('biom')] = 0
ksiin[metabolites.index('putr')] = 0
ksiin[metabolites.index('3HB')] = 0

print('Flux entrant ', 'biom', ' :', ksiin[metabolites.index('biom')])
print('Flux entrant ', 'putr', ' :', ksiin[metabolites.index('putr')])
print('Flux entrant ', '3HB', ' :', ksiin[metabolites.index('3HB')])
print()


ksi0[metabolites.index('biom')] = biom_init * convmg_L__CFU_gFW

print('Initial biomass (mg/L)', ksi0[metabolites.index('biom')])



#time for the ODEs
t = np.linspace(t_start*24,  t_end*24, t_step)

# solve the ODEs
y = odeint(MM, ksi0, t)  # integrate

# remove negative values for other substrates
for (ind_othersubstrate, other_substrate) in enumerate(other_substrates):
    for ind_time, time in enumerate(t):
        if y[ind_time, nb_substrate + ind_othersubstrate] < 0:
            # print(other_substrate, " ", time, " ", y[ind_time, nb_substrate + ind_othersubstrate])
            y[ind_time, nb_substrate + ind_othersubstrate] = 0



#print biomass at 72 hours
ind_time_72 = 0

for ind_time, time in enumerate(t):
    if time >= 72:
        ind_time_72 = ind_time
        break

print('Biomass at 72 hours: ', y[ind_time_72, metabolites.index('biom')])



#########################
# #compute error
#######################

errors = np.zeros(1)

for (ind_metabolite, metabolite) in enumerate(metabolites):

    error = np.zeros(1)

    if np.nanmean(data_exp[metabolite]) != 0.0:
        #print(metabolite)

        # compute error compared to experimental data
        y_time_metabolite = []

        for datatime in data_exp['Time']:
            for ind_time, time in enumerate(t):
                if time >= datatime * 24:
                    #print(time)
                    if metabolite == 'biom':
                        try:
                            y_time_metabolite = y_time_metabolite + [math.log10(y[ind_time, ind_metabolite]/ convmg_L__CFU_gFW)]
                        except:
                            print("problem", y[ind_time, ind_metabolite])
                            y_time_metabolite = 0
                    elif metabolite == '3HB':
                        y_time_metabolite = y_time_metabolite + [y[ind_time, ind_metabolite]/104.105 * nb_C[ind_metabolite]]
                    else:
                        y_time_metabolite = y_time_metabolite + [y[ind_time, ind_metabolite] * nb_C[ind_metabolite]]
                    break

        error = np.concatenate((error, 1/np.nanmean(data_exp[metabolite]) * (y_time_metabolite - data_exp[metabolite])))

    errors = np.concatenate((errors, error))

squared_error = np.nansum(np.square(errors))

print()
print(f'Current global squared error:', squared_error)




#########################
# #plot results
#######################
tot_plots = len(metabolites) + 1

fig, axs = plt.subplots(math.ceil(tot_plots/2), 2, figsize=(25, 20))

count_plot = 0


for (ind_metabolite, metabolite) in enumerate(metabolites):

    if ksi0[ind_metabolite] != 0.0 or metabolite == 'putr' or metabolite == '3HB':

        # plot results

        ax = axs[count_plot%math.ceil(tot_plots/2), count_plot//math.ceil(tot_plots/2)]
        if metabolite == 'biom':

            ax.plot(t / 24, y[:, ind_metabolite]/1000, 'b-')
            ax.set_ylabel('Biom (g/L)')
            ax.errorbar(data_exp['Time'], [pow(10, value)*convmg_L__CFU_gFW/1000 for value in data_exp[metabolite]], fmt='bo',
                            yerr=[pow(10, value)*convmg_L__CFU_gFW/1000 for value in data_exp_std[metabolite]])
            ax.set_xlabel('time (d)')
            ax.set_title(metabolite)
            ax.grid(True)
            ax.set_ybound(lower=0)

            count_plot = count_plot + 1

            ax = axs[count_plot % math.ceil(tot_plots / 2), count_plot // math.ceil(tot_plots / 2)]
            ax.plot(t / 24, [math.log10(value / convmg_L__CFU_gFW) for value in y[:, ind_metabolite]], 'b-')
            ax.set_ylabel('cFU/gFW')
            ax.errorbar(data_exp['Time'], data_exp[metabolite], fmt='bo',
                        yerr=data_exp_std[metabolite])


        elif metabolite == '3HB':
            ax.plot(t / 24, y[:, ind_metabolite]/104.105, 'b-')
            ax.set_ylabel(metabolite + ' (mM)')
            ax.errorbar(data_exp['Time'], data_exp[metabolite] / nb_C[ind_metabolite], fmt='bo',
                        yerr=data_exp_std[metabolite] / nb_C[ind_metabolite])

        else:
            ax.plot(t/24, y[:, ind_metabolite], 'b-')
            ax.set_ylabel(metabolite + ' (mM)')
            ax.errorbar(data_exp['Time'], data_exp[metabolite]/nb_C[ind_metabolite], fmt='bo', yerr=data_exp_std[metabolite]/nb_C[ind_metabolite])


        ax.set_xlabel('time (d)')
        ax.set_title(metabolite)
        ax.grid(True)
        ax.set_ybound(lower=0)

        count_plot = count_plot + 1



plt.tight_layout()
plt.savefig('output/png/SimuXylemContinue.png')
#plt.show()

plt.clf()



#plot interpolated data

MEDIUM_SIZE = 20

plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels

fig, axs = plt.subplots(4, 1, figsize=(12, 14))

ax = axs[0]
ax.plot(ts, transpis, 'b-')
ax.errorbar(data_transpi['dpi'], data_transpi['transpi_tot (mL)'], fmt='bo', yerr=data_transpi['std transpi_tot (mL)'])
ax.set_ylabel('Transpiration (mL/d)')
ax.set_xlabel('Time (d)')
#ax.set_title('Transpiration (mL/d)')
ax.grid(True)

ax = axs[1]
ax.plot(ts, fws, 'b-')
ax.errorbar(data_fw['Time (dpi)'], data_fw['fwplante (g)'], fmt='bo', yerr=data_fw['std_fwplante (g)'])
ax.set_ylabel('Plant fresh weight (g)')
ax.set_xlabel('Time (d)')
#ax.set_title('Fresh weight (g)')
ax.grid(True)

ax = axs[2]
ax.plot(ts, vxyls, 'b-')
ax.plot(data_fw['Time (dpi)'], [percentage_xyl * fw_plante / densite_plante for fw_plante in data_fw['fwplante (g)']], 'bo')
ax.set_ylabel('Volume Xylem (mL)')
ax.set_xlabel('Time (d)')
#ax.set_title('Volume Xylem (mL)')
ax.grid(True)

ax = axs[3]
ax.plot(ts, Ds, 'b-')
ax.set_ylabel('Dilution rate (h-1)')
ax.plot(data_fw['Time (dpi)'], [data_fw['transpi_tot'][i]/24 * 1 / (percentage_xyl * fw_plante / densite_plante) for i, fw_plante in enumerate(data_fw['fwplante (g)'])], 'bo')
ax.set_xlabel('Time (d)')
#ax.set_title('Dilution rate (h-1)')
ax.grid(True)

ax.set_ybound(lower=0)

plt.tight_layout()

plt.savefig('output/png/SimuXylemContinueInterpData.png')




#########################
# #Export to excel
#######################
df = pd.DataFrame({'time (d)': t/24})

for (ind_metabolite, metabolite) in enumerate(metabolites):
    if metabolite =='3HB':
        df['3HB (mM)'] = y[:, ind_metabolite] / 104.105
    elif metabolite == 'biom':
        df['Biom (g/L)'] = y[:, ind_metabolite] / 1000
        df['cFU/gFW'] = [math.log10(value / convmg_L__CFU_gFW) for value in y[:, ind_metabolite]]
    else:
        df[metabolite + ' (mM)'] = y[:, ind_metabolite]


df.to_excel('output/xls/XylemContinu_metabolites.xlsx',
             sheet_name='XylemContinu_metabolites', index=False)

