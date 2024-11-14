import numpy as np
import matplotlib.pyplot as plt
from wavelen_work import *
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import os
import re


def get_nearest(x_query, x_vals):
    nearest_index = np.argmin(np.abs(x_vals - x_query))
    return nearest_index


def chi_squared(synth, observed):
    return np.sum((synth - observed) ** 2 / observed)

def chi_2(path, synthetics, stardata, k, spectral_lines,name=None, start=None, end=None, plot=None):
    if k < 18500:
        j="h"
    else :
        j="k"
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + j), stardata.get("v_" + j)), stardata.get("flux_" + j), k)

    m = len(normal['flux_normalised']) // 2

    # Longueurs d'onde observées
    observed_wavelengths = np.array(normal['z_wavelen'][m - normal['index_closest_with_tolerance']: m + normal['index_closest_with_tolerance']+1])
    # Flux observé : portion autour de l'indice central
    observed = np.array(normal['flux_normalised'][m - normal['index_closest_with_tolerance']: m + normal['index_closest_with_tolerance']+1])
    
    if start is not None and end is not None:
        observed_wavelengths = np.array(normal['z_wavelen'][get_nearest(start,np.array(normal['z_wavelen'])):get_nearest(end,np.array(normal['z_wavelen']))])
        observed = np.array(normal['flux_normalised'][get_nearest(start,np.array(normal['z_wavelen'])):get_nearest(end,np.array(normal['z_wavelen']))])
    

        
    taille = (observed_wavelengths[-1]-observed_wavelengths[0])/2
    # Création de la liste pour stocker les spectres 
    synthetic_spectra = []
    observed_spectra =[]
    observed_spectra_w =[]
    synthetic_spectra_w = []

    # Calcul des spectres synthétiques
    for syntha in synthetics:
            
            synt_flux = zoom_syntspec(path, syntha, k, taille)["synt_flux"]
            synt_flux_w = zoom_syntspec(path, syntha, k, taille)["synt_wavelen"]
            synthetic_spectra.append(synt_flux)
            synthetic_spectra_w.append(synt_flux_w)

            interp_func = interp1d(observed_wavelengths, observed, kind='cubic')
            xre = np.linspace(observed_wavelengths[0], observed_wavelengths[-1], len(synt_flux))
            yre = interp_func(xre)

            observed_spectra.append(yre)
            observed_spectra_w.append(xre)

    # Conversion en tableau numpy
    observed_spectra = np.array(observed_spectra)
    synthetic_spectra = np.array(synthetic_spectra)
    # results = np.array(results)

    #Calcul des valeurs de chi carré
    chi_squared_values = []
    for i in range(len(observed_spectra)): 
        chi_squared_values.append(chi_squared(synthetic_spectra[i], observed_spectra[i]))

    items = list(synthetics.items())

    #Affichage des résultats
    # for i, chi_squared_value in enumerate(chi_squared_values):
    #     element = items[i]
    #     print(f"Chi-carré pour le spectre synthétique de la raie {k} pour {element[1]} : {chi_squared_value}")

#tracé du graphe
    if plot is True:
        f = plt.figure(figsize=(15,8), dpi = 400)
        gs = f.add_gridspec(2, hspace=0.2)
        ax = gs.subplots(sharex=False, sharey=True)

        ax[0].set_xlim(k-11,k+11)
        ax[1].set_xlim(k-3,k+3)

        ax[0].xaxis.set_major_locator(MultipleLocator(5))
        ax[0].xaxis.set_minor_locator(MultipleLocator(1))
        ax[1].xaxis.set_major_locator(MultipleLocator(1))
        ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))

        # Récupération de toutes les longueurs d'onde dans une liste unique avec les éléments associés
        all_wavelengths = []
        for element, wavelengths in spectral_lines.items():
            for wavelength in wavelengths:
                all_wavelengths.append((wavelength, element))

        # Tri de la liste des longueurs d'onde
        all_wavelengths.sort()  # Trie par longueur d'onde en ordre croissant

        # Boucle pour tracer les lignes et afficher les éléments
        text_height_base = 1.2
        text_height = text_height_base  # Initialisation de la hauteur actuelle

        for i, (z, element) in enumerate(all_wavelengths):
            # Détection de la proximité avec la longueur d'onde précédente
            if i > 0 and abs(z - all_wavelengths[i - 1][0]) < 0.7:
                text_height += 0.07  # Si trop proche du précédent, augmenter la hauteur de 0.05
            else:
                text_height = text_height_base  # Réinitialisation de la hauteur si suffisamment éloigné

            # Vérification des conditions pour tracer les lignes
            if k - 11 <= z <= k + 11:
                # Tracé de la ligne verticale
                ax[0].axvline(x=z, ymin=0.7, ymax=0.8, color='black', linewidth=0.5)
                # Ajout du texte de l'élément avec hauteur adaptée
                ax[0].text(z, text_height, s=element, color='black', fontsize=10, ha='center')
            
            if k - 3 <= z <= k + 3:
                # Tracé de la ligne verticale
                ax[1].axvline(x=z, ymin=0.7, ymax=0.8, color='black', linewidth=0.5)
                # Ajout du texte de l'élément avec hauteur adaptée
                ax[1].text(z, text_height, s=element, color='black', fontsize=10, ha='center')
        for ax in ax :
            ax.plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=1, color='black', label="Spectre observé : " + stardata.get("starname"))
            
            if name is not None:
                ax.axvline(x=k, ymin=0.7, ymax=0.8, color='black', linewidth=0.5)
                ax.text(k, 1.2, s=name, color='black', fontsize=10, ha='center')
            
            for synth in synthetics : 
                AX2 = syntspec(path+synth)
                ax.plot(AX2['wavelen'], AX2['flux'], linewidth = 1, label="Synth : "+synthetics.get(synth))

            # ax.plot(observed_spectra_w[1],observed_spectra[1])

            ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
            ax.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
            ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
            ax.set_ylim(0., 1.4)
            ax.tick_params(axis = 'both', labelsize = 10)
            ax.legend(loc = "lower left", fontsize = 10)
            ax.axvspan(observed_wavelengths[0], observed_wavelengths[-1], color='bisque', alpha=0.3)
            ax.set_xlabel("Longueur d'onde (Å)", fontsize  = 10)
            ax.set_ylabel("Flux normalisé", fontsize  = 10)
        plt.show()

    return {"chi_squared_values":chi_squared_values}


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

def chi_plot_ABU(abu, chi_squared,element, k, raie, raie_OH):
    plt.switch_backend('pgf')
    abu = np.array(abu)
    # Perform quadratic fitting
    params, _ = curve_fit(quadratic, abu, chi_squared)
    a, b, c = params  # coefficients of the quadratic fit
    fit_x = np.linspace(abu.min(), abu.max(), 100)
    fit_y = quadratic(fit_x, a, b, c)

    # Find minimum of the quadratic fit
    min_log_e = -b / (2 * a)  # Vertex of the parabola
    min_chi_squared = quadratic(min_log_e, a, b, c)

    # Plot
    f = plt.figure(figsize=(8, 6), dpi = 400)
    gs = f.add_gridspec(1, hspace=0.2)
    ax = gs.subplots(sharex=False, sharey=True)
    ax.scatter(abu, chi_squared, color="darkblue", marker="x", label=f"Raie de {raie} en {k} Å")  # Data points
    # ax.plot(fit_x, fit_y, color="lightgray", label="Ajustement quadratique")  # Quadratic fit line
    ax.scatter(min_log_e, min_chi_squared, color="red", marker="x", label=f"Minimum en {min_log_e:.2f}")  # Minimum point
    ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
    ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
    ax.tick_params(axis = 'both', labelsize = 16)
    # Labeling
    ax.set_xlabel(f"$\\log \\epsilon_{{\\mathrm{{{element}}}}}$", fontsize=16)
    ax.set_ylabel(r"$\chi^2$", fontsize=16)
    plt.legend(loc="upper left", fontsize=16)
    plt.show()
    plt.savefig('graphique_test.pgf')
    raie_OH[k] = min_log_e

def interpolated_chi_squared(x, abundances, macroturbulences, chi_squared_values):
    ab, mt = x
    return griddata((abundances, macroturbulences), chi_squared_values, (ab, mt), method='cubic')

# Fonction de minimisation et visualisation
def double_chi(abundances, macroturbulences, chi_squared_values):
    ab_min, ab_max = min(abundances), max(abundances)
    mt_min, mt_max = min(macroturbulences), max(macroturbulences)

    # Générer une grille fine pour la visualisation
    grid_ab, grid_mt = np.mgrid[ab_min:ab_max:100j, mt_min:mt_max:100j]
    grid_chi_squared = griddata((abundances, macroturbulences), chi_squared_values, (grid_ab, grid_mt), method='cubic')

    # Estimation initiale pour la minimisation
    initial_guess = [np.mean(abundances), np.mean(macroturbulences)]
    
    # Minimisation avec interpolation du chi carré
    result = minimize(
        interpolated_chi_squared, initial_guess, 
        args=(abundances, macroturbulences, chi_squared_values),  # Arguments supplémentaires
        bounds=[(ab_min, ab_max), (mt_min, mt_max)]
    )

    # Récupérer les valeurs minimales estimées
    best_abundance = result.x[0]
    best_macroturbulence = result.x[1]
    min_chi_squared_value = result.fun

    print("Meilleure abondance estimée :", best_abundance)
    print("Meilleure macroturbulence estimée :", best_macroturbulence)
    print("Valeur minimale de χ² estimée :", min_chi_squared_value)

    # Visualisation
    plt.figure(figsize=(10, 8))
    plt.contourf(grid_ab, grid_mt, grid_chi_squared, levels=20, cmap='magma')
    plt.colorbar(label='$\chi^2$')
    plt.scatter(abundances, macroturbulences, c=chi_squared_values, cmap='magma', edgecolor='orange', label='Données')
    plt.plot(best_abundance, best_macroturbulence, 'rx', markersize=10, label='Min $\chi^2$ interpolé')
    plt.xlabel("Abondance")
    plt.ylabel("Macroturbulence")
    plt.legend()
    plt.show()


def plot_chi_squared_MAC(path, stardata,lines_BD22, raie_propre, dossier_element):
    abondance = []
    macroturbulence = []
    synth = {}
    target_path = path+dossier_element+"/"
    for filename in os.listdir(target_path):
        # print(filename)
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":

            match = re.search(r'MAC_-(.*?)_', filename)
            # print(match)
            match2 = re.search(r'abu_(.*?).conv', filename)
            
            if match and match2:
                chars_after_MAC = match.group(1)
                # print(chars_after_MAC)
                chars_after_ABU = match2.group(1)

            synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU} & $v_{{\\mathrm{{macro}}}}$ = {chars_after_MAC} $km s^{{-1}}$"
            abondance.append(float(chars_after_ABU))
            macroturbulence.append(float(chars_after_MAC))

    # Calculer les valeurs de chi_squared pour chaque raie
    chi_squared_values = None  # Initialiser au cas où aucune valeur ne soit calculée
    for raie in raie_propre:
        for k in raie_propre[raie]:
            result = chi_2(target_path, synth, stardata, k, lines_BD22, raie)
            chi_squared_values = result.get("chi_squared_values", [])  # Assurez-vous de récupérer chi_squared_values

    # Appeler double_chi si chi_squared_values a bien été défini
    if chi_squared_values is not None:
        double_chi(abondance, macroturbulence, chi_squared_values)
    else:
        print("Erreur: chi_squared_values n'a pas été calculé.")



def plot_chi2_simple_MAC(path, dossier_element, stardata, raie_propre, lines_BD22, target_pairs=None):
    """
    Affiche les graphes pour les spectres synthétiques qui correspondent aux paires d'abondance et de macroturbulence cibles.

    target_pairs : Liste de tuples avec des paires d'(abondance, macroturbulence). Exemple : [(6.5, 10.5), (7.0, 9.0)]
    """
    synth = {}
    target_path = path+dossier_element+"/"
    
    for filename in os.listdir(target_path):
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":
            # Rechercher les valeurs de MAC et d'abondance dans le nom du fichier
            match = re.search(r'MAC_-(.*?)_', filename)
            
            match2 = re.search(r'abu_(.*?).conv', filename)
            
            if match and match2:
                chars_after_MAC = match.group(1)
                chars_after_ABU = match2.group(1)

                # Vérifier si la paire (abondance, macroturbulence) est dans target_pairs
                if target_pairs is None or (float(chars_after_ABU), float(chars_after_MAC)) in target_pairs:
                    synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU} & $v_{{\\mathrm{{macro}}}}$ = {chars_after_MAC} $km s^{{-1}}$"
                    
    # Parcourir les raies et appeler la fonction chi_2 pour chaque raie spécifique
    for raie in raie_propre:
        for k in raie_propre[raie]:
            chi_2(target_path, synth, stardata, k, lines_BD22, plot=True)


def plot_chi2_simple_ABU(path, dossier_element, stardata, raie_propre, lines_BD22,start=None,end=None, target_pairs=None):
    """

    """
    synth = {}
    target_path = path+dossier_element+"/"
    
    for filename in os.listdir(target_path):
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":
            # Rechercher les valeurs de MAC et d'abondance dans le nom du fichier
            
            match = re.search(r'abu_(.*?).conv', filename)
            
            if match:
                chars_after_ABU = match.group(1)

                # Vérifier si la paire (abondance, macroturbulence) est dans target_pairs
                if target_pairs is None or (float(chars_after_ABU)) in target_pairs:
                    synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU}"
                    
    # Parcourir les raies et appeler la fonction chi_2 pour chaque raie spécifique
    for raie in raie_propre:
        for k in raie_propre[raie]:
            chi_2(target_path, synth, stardata, k, lines_BD22,start,end, plot=True)