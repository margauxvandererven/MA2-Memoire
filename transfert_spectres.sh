#!/bin/bash

# Définir l'adresse du serveur et les répertoires
serveur="student@10.149.16.137"
script_com="/home/student//SPECTRUM/Turbospectrum_NLTE-master/COM/script-NLTE-IR_BD-221742.com"
dossier_dist="/home/student//SPECTRUM/Turbospectrum_NLTE-master/COM/syntspec/"
dossier_local="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/CO_mol/"


# Intervalle de temps pour les fichiers récents (par exemple, 1 heure)
intervalle="1"  # heures

# Connexion SSH et exécution du script .com
echo "Exécution du script sur le serveur distant..."
ssh "$serveur" "cd $(dirname "$script_com"); ./$(basename "$script_com")"

# Rechercher et télécharger les fichiers générés au cours de la dernière heure
echo "Téléchargement des fichiers générés récemment..."

# Utiliser SSH pour trouver les fichiers modifiés récemment et les transférer via SCP
fichiers_recents=$(ssh "$serveur" "find $dossier_dist -type f -name '*.conv' -mmin -$((intervalle * 30))")

# Boucle pour télécharger chaque fichier récent
for fichier in $fichiers_recents; do
    scp "$serveur:$fichier" "$dossier_local"
done

# Message de confirmation
echo "Les spectres récents ont été téléchargés avec succès."
