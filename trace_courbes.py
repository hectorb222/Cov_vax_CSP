import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def importation(filename):
    """ La fonction convertit les donnes csv en un tableau de valeurs type DataFrame
    Parameters
    ----------
    filename: str, default='./data/fr-esr-parcoursup.csv'
        nom du fichier csv avec les donnees a traiter

    Returns
    -------
    tableau: DataFrame
        tableau de valeurs du fichier importe en csv convertit en DataFrame
    """

    return pd.read_csv(filename, sep=",")


array_dead = np.zeros((6, 5))
array_infected = np.zeros((6, 5))

# array_dead[i][j] := nombre de morts cumulés pour strat i et scénar j
# array_infected[i][j] := nombre de morts cumulés pour strat i et scénar j

df_1 = importation("./dataframe3.csv")
res = []

for t in range(0, 3540, 118):
    res.append(df_1.iloc[:, t:t+118])

# Numero de la stratégie x scénario
for i in range(6):
    for j in range(5):
        n = 5*i + j

        df_1_sains = res[n].iloc[:, 1: 28]
        df_1_infected = res[n].iloc[:, 55: 82]
        df_1_retablis = res[n].iloc[:, 82: 109]
        df_1_mort = res[n].iloc[:, 109: 136]

        list_infected = df_1_infected.sum(axis=1)
        list_infected_cum = list_infected.cumsum()


# Courbe total des rétablis
        list_retablis = df_1_retablis.sum(axis=1)

# Courbe total des sains
        list_sains = df_1_sains.sum(axis=1)

# Deces cumulé
        list_mort = df_1_mort.sum(axis=1)
        array_dead[i][j] = list_mort[364]
        array_infected[i][j] = list_mort[364] + list_retablis[364]


array_dead = np.matrix.transpose(array_dead)
array_infected = np.matrix.transpose(array_infected)

array_infected_vf = np.zeros((4, 6))
array_dead_vf = np.zeros((4, 6))

for i in range(0, 3):
    array_dead_vf[i] = array_dead[i]
    array_infected_vf[i] = array_infected[i]
array_dead_vf[3] = array_dead[4]
array_infected_vf[3] = array_infected[3]


fig, ax = plt.subplots(1, 1)

# img = ax.imshow(array_dead_vf, cmap='YlGnBu', extent=[-1, 1, -1, 1])
# plt.title("Nombre de décès par couple stratégie/scénario")

img = ax.imshow(array_infected_vf, cmap='BuGn', extent=[-1, 1, -1, 1])
plt.title("Nombre d'infections par couple stratégie/scénario")

x_label_list = ['Témoin', 'Plus vulnérables', 'Acteurs économiques',
                'Travailleurs essentiels', 'Non télétravailleurs', 'Très exposés']
y_label_list = ['Aucune restriction', 'Déconfinement progressif',
                'Confinement partiel', 'Confinement total']


ax.set_xticks(np.linspace(-1+1/6, 1-1/6, 6))
ax.set_yticks(np.linspace(-1+1/4, 1-1/4, 4))

ax.set_xticklabels(x_label_list)
ax.set_yticklabels(y_label_list)

plt.xticks(rotation=30, ha="right")
# plt.yticks(rotation=45, ha="right")

# plt.xlabel("Prioritaires à la vaccination")
# plt.ylabel("Scénario")

plt.colorbar(img)


plt.show()
