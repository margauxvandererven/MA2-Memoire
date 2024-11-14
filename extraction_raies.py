import numpy as np
from wavelen_work import *
from readmultispec import *
import json

starname = input("Entrez le nom de l'étoile (ex. BD-221742) : ")

# lecture des spectres fichiers .fits
read_spec_h = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/H_Sneden/"+starname+"H.cont.fits")
read_spec_k = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/K_Sneden/"+starname+"K.cont.fits")

stardata = {
    "starname"  : str(starname),
    "wavelen_h" : list(read_spec_h['wavelen']),
    "flux_h" : read_spec_h['flux'], 
    "wavelen_k" : list(read_spec_k['wavelen']), 
    "flux_k" : read_spec_k['flux'], 
    "v_h" : 176.8 * 1e3, # redshift déterminé pour la bande H
    "v_k" : 176.8* 1e3 # redshift déterminé pour la bande K
}


def get_nearest_y(x_query, x_vals, y_vals):
    nearest_index = np.argmin(np.abs(x_vals - x_query))
    return y_vals[nearest_index]


def extraction_raies(pathtofichiertxt,pathsynth, stardata):
    """
    Extrait les raies qui sont prononcées à 5% dans le spectre observé et atomique
    """
    with open(pathtofichiertxt,"r",encoding="utf8") as file :
        lines = file.readlines()

    lines_star = {
    }
    W_lambda = {}

    for line in lines : 
        split = line.split()
        element = str(split[5][1:])
        lines_star[element.capitalize()+" "+ split[6]] = []
        W_lambda[element.capitalize()+" "+ split[6]] = []

    for line in lines : 
        split = line.split()
        element = split[5][1:] 
        if element in split[5]:
            raies = lines_star.get(element.capitalize()+" "+ split[6])
            largeur = W_lambda.get(element.capitalize()+" "+ split[6])

            if 14600 < float(split[0]) < 18400 : 
                normal = normalisation(redshift_wavelen(stardata.get("wavelen_h"), stardata.get("v_h")),stardata.get("flux_h"),float(split[0]))
                AX = syntspec(pathsynth + "atomic")
                y_nearest = get_nearest_y(float(split[0]), np.array(normal['z_wavelen']), np.array(normal['flux_normalised']))
                y_nearest_synth = get_nearest_y(float(split[0]), np.array(AX['wavelen']), np.array(AX['flux']))
                if y_nearest < 0.95 and y_nearest_synth < 0.95 : 
                    raies.append(float(split[0]))
                    largeur.append(float(split[3])*10**8)

            if 18700 < float(split[0]) < 25300 : 
                normal = normalisation(redshift_wavelen(stardata.get("wavelen_k"), stardata.get("v_k")),stardata.get("flux_k"),float(split[0]))
                AX = syntspec(pathsynth + "atomic2")
                y_nearest = get_nearest_y(float(split[0]), np.array(normal['z_wavelen']), np.array(normal['flux_normalised']))
                y_nearest_synth = get_nearest_y(float(split[0]), np.array(AX['wavelen']), np.array(AX['flux']))
                if y_nearest < 0.95 and y_nearest_synth < 0.95 : 
                    raies.append(float(split[0]))
                    largeur.append(float(split[3])*10**8)
            else:
                continue

    return {"lines_star":lines_star, "W_lambda":W_lambda}

path="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/"+starname+"b/"
lines_BD22 = extraction_raies(path+"raie_atom.txt",
                              path, stardata)["lines_star"]
W_lambda = extraction_raies(path+"raie_atom.txt",
                              path, stardata)["W_lambda"]

def filter_close_wavelengths(wavelengths, threshold=0.1):
    """
    Lorsque plusieurs raies du même élément se collent, on ne garde qu'une raie à afficher
    """
    if not wavelengths:
        return []
    # Trier les longueurs d'onde
    wavelengths = sorted(wavelengths)
    filtered_wavelengths = []
    group = [wavelengths[0]]
    
    for i in range(1, len(wavelengths)):
        # Vérifier si la différence avec la précédente est inférieure au seuil
        if wavelengths[i] - wavelengths[i - 1] < threshold:
            group.append(wavelengths[i])
        else:
            # Ajouter la valeur médiane du groupe à la liste finale
            filtered_wavelengths.append(np.median(group))
            group = [wavelengths[i]]
    
    # Ajouter le dernier groupe si présent
    if group:
        filtered_wavelengths.append(np.median(group))
    
    return filtered_wavelengths

# Application du filtre sur chaque liste d'élément dans lines_BD22
for element, wavelengths in lines_BD22.items():
    lines_BD22[element] = filter_close_wavelengths(wavelengths)

with open("extraction_raies_"+stardata.get("starname")+".txt", "w") as fichier:
    json.dump(lines_BD22, fichier, indent=4, ensure_ascii=False)

# Dictionnaire de longueurs d'onde filtrées
# Assurez-vous que les longueurs d'onde filtrées sont dans le même ordre que celles d'origine
filtered_equivalent_widths = {}

# Itérer sur chaque élément pour filtrer les largeurs équivalentes
for element, filtered_wavelengths in lines_BD22.items():
    original_wavelengths = lines_BD22[element]  # Liste des longueurs d'onde d'origine
    original_widths = W_lambda[element]  # Liste des largeurs équivalentes d'origine

    # Construire une liste filtrée pour les largeurs équivalentes correspondant aux longueurs d'onde
    filtered_widths = [
        original_widths[original_wavelengths.index(wav)] 
        for wav in filtered_wavelengths if wav in original_wavelengths
    ]

    # Ajouter à notre dictionnaire final
    filtered_equivalent_widths[element] = filtered_widths


# Étape 1 : Combiner les longueurs d'onde et les largeurs équivalentes dans une seule liste
combined_list = []
for element, wavelengths in lines_BD22.items():
    widths = filtered_equivalent_widths[element]
    combined_list.extend([(wavelength, width, element) for wavelength, width in zip(wavelengths, widths)])

# Étape 2 : Trier la liste par longueur d'onde
combined_list.sort()

# Étape 3 : Parcourir la liste triée et appliquer la règle de comparaison
filtered_combined_list = []
i = 0

while i < len(combined_list):
    current_wavelength, current_width, current_element = combined_list[i]
    # Vérifier si le prochain élément est à moins de 0.5 de distance
    if i < len(combined_list) - 1 and abs(combined_list[i + 1][0] - current_wavelength) < 0.5:
        # Comparer avec la longueur d'onde suivante
        next_wavelength, next_width, next_element = combined_list[i + 1]
        # Garder celle avec la largeur équivalente la plus grande
        if current_width >= next_width:
            filtered_combined_list.append((current_wavelength, current_width, current_element))
        else:
            filtered_combined_list.append((next_wavelength, next_width, next_element))
        # Passer au-delà de la paire comparée
        i += 2
    else:
        # Si aucune proximité détectée, conserver l'élément courant
        filtered_combined_list.append((current_wavelength, current_width, current_element))
        i += 1

# Étape 4 : Reconstruire le dictionnaire final
final_filtered_lines_BD22 = {}
final_filtered_equivalent_widths = {}

for wavelength, width, element in filtered_combined_list:
    if element not in final_filtered_lines_BD22:
        final_filtered_lines_BD22[element] = []
        final_filtered_equivalent_widths[element] = []
    final_filtered_lines_BD22[element].append(wavelength)
    final_filtered_equivalent_widths[element].append(width)