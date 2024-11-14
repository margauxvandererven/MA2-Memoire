import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from wavelen_work import *

# lines = {
#     "S": [15478.48, 22507.60],
#     "Na": [16373.87 , 16388.85 , 22056.43 , 22083.66  , 23379.14],
#     "Al": [16718.97 , 16750.60 , 16763.37 , 17699.05 , 21093.08 , 21163.80 , 21208.18], 
#     "Zn": [16505.18],
#     "Si": [16434.93, 20804.20 , 20890.37  , 20926.14],
#     "K": [15163.09 , 15168.40],
#     "Ca": [16150.76 , 16155.24 , 16157.36 ,  20962.57 , 20972.53 , 21113.90 , 22607.94 , 22624.96 , 22626.72 , 22651.18 , 22653.58],
#     "Sc": [21730.36 , 21812.24 , 21842.70],
#     "Ti": [16330.54 , 21149.62],
#     "V": [15924.81],
#     "Cr": [15680.06 , 15860.21 , 17708.73 ], 
#     "Mn": [15217.74 , 15262.49],
#     "Co": [16757.64],
#     "Ni": [16310.50 , 16363.09 , 16815.46 , 16818.74 , 16867.28 , 17306.52 , 20957.14 , 21167.93 , 21570.06 , 21945.50 , 22596.93 ],
#     "Cu": [16005.64, 16638.98], 
#     "Zn": [16505.18],
#     "Y": [21260.45 , 22543.84],
#     "Ce": [15277.65 , 15829.83 , 15977.12 , 16595.18],
#     "Nd": [15368.14 , 16053.63 , 16262.04],
#     "Yb": [16498.40],
#     "F": [22699.49 , 22714.59 ,  22778.25 , 22826.86 , 22886.73 , 22957.94],
#     "Mg" : [21059.76, 21060.89, 21458.87]
# }



# fonction permettant de tracer plusieurs plot à la suite suivant l'analyse demandée

def plot_lines(ax, l, path, synthetics, stardata, n, k, i, m, j, taille, spectral_lines):
    """
    Ouvre et trace le spectre observé normalisé. Trace les spectres synthétiques "synthetics". Le tout centré sur les raies séléctionnées
    avec une largeur demandée (taille)
    ax : le plot
    l : ?
    path : chemin des spectres synthétiques
    synthetics : spectre synthétiques
    stardata : données de l'étoile
    n : nombre d'éléments dans la liste de raies
    k : raie de centrage
    i : nom de la raie
    m : bande h ou k
    j : "", "2" ou "3"
    taille : taille du zoom
    """
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)),stardata.get("flux_" + m),k)
    if n == 1:
        ax.plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='black', label="observé : " + stardata.get("starname"))
        for mol in synthetics:
            AXb = syntspec(path + mol + j)
            ax.plot(AXb['wavelen'], AXb['flux'], label= "Synth: " + str(synthetics.get(mol)), linewidth=0.8)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_tick_params(direction='in', length=10, which='major')
        ax.xaxis.set_tick_params(direction='in', length=6, which='minor')
        for element, wavelengths in spectral_lines.items():
            for z in wavelengths:
                if k - taille <= z <=  k + taille:
                    # Tracé de la ligne verticale
                    ax.axvline(x=z, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
                    # Ajout du texte de l'élément au-dessus de la ligne
                    ax.text(z, 1.2, s=element, color='black', fontsize=12, ha='center')
        ax.text(k, 1.2, s=i, color='black', fontsize=12, ha='center')
        ax.axvline(x=k, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
        ax.set_xlim(k - taille, k + taille)
        ax.set_xlim(k - taille, k + taille)
        ax.set_ylim(0.2, 1.3)

        ax.legend(loc='lower right')
        plt.savefig('mon_graphique.pdf', format='pdf')
    else:
        ax[l].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='black', label="observé : " + stardata.get("starname"))
        for mol in synthetics:
            AXb = syntspec(path + mol + j)
            ax[l].plot(AXb['wavelen'], AXb['flux'], label="Synth: " + str(synthetics.get(mol)), linewidth=0.8)
        ax[l].xaxis.set_major_locator(MultipleLocator(5))
        ax[l].xaxis.set_minor_locator(MultipleLocator(1))
        ax[l].xaxis.set_tick_params(direction='in', length=10, which='major')
        ax[l].xaxis.set_tick_params(direction='in', length=6, which='minor')
        for element, wavelengths in spectral_lines.items():
            for z in wavelengths:
                if k - taille <= z <=  k + taille:
                    # Tracé de la ligne verticale
                    ax[l].axvline(x=z, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
                    # Ajout du texte de l'élément au-dessus de la ligne
                    ax[l].text(z, 1.2, s=element, color='black', fontsize=12, ha='center')
        ax[l].axvline(x=k, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
                    # Ajout du texte de l'élément au-dessus de la ligne
        ax[l].text(k, 1.2, s=i, color='black', fontsize=12, ha='center')
        ax[l].set_xlim(k - taille, k + taille)
        ax[l].set_ylim(0.2, 1.3)
        ax[l].legend(loc='lower right')
        plt.savefig('mon_graphique.pdf', format='pdf')


def zoom_lines(lines, path, synthetics, stardata, taille_zoom, spectral_lines) :
    """
    lines : raies que l'on souhaite observer
    path : chemin vers spectres synthétiques
    synthetics : spectres synthétiqyes
    stardata : données de l'étoile
    taille_zoom : taille du zoom sur chaque raie
    """
    for i in list(lines.keys()):
        n = len(lines.get(i))
        y_size = n*3
        f = plt.figure(figsize=(15,y_size), dpi = 200)
        gs = f.add_gridspec(n, hspace=0.2)
        ax = gs.subplots(sharex=False, sharey=True)
        for k in lines.get(i):
            l = lines.get(i).index(k)
            if  k < 18500 : 
                plot_lines(ax,l, path, synthetics, stardata, n, k, i, "h", "", taille_zoom, spectral_lines)
            else:
                plot_lines(ax,l, path, synthetics, stardata, n, k, i, "k", "2", taille_zoom, spectral_lines)

            # plt.savefig("/Users/margauxvandererven/Documents/unif2023-2024/spectre_IR/output/"+stardata.get("starname")+"/"+i)


def zoom_lines_analyse(ax, path, synthetics, molecules, stardata, j, m, k, i) : 
    # AX = syntspec(path+filename)

    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)),stardata.get("flux_" + m),k)

    ax[0].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=1, color='black', label="Spectre observé : " + stardata.get("starname"))
    ax[1].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=1, color='black', label="Spectre observé : " + stardata.get("starname"))
    # ax[2].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='lightgray', label="Spectre observé : " + stardata.get("starname"))

    for mol in molecules:
        AXb = syntspec(path + mol + j)
        ax[0].plot(AXb['wavelen'], AXb['flux'], label="Synth : atom +" + molecules.get(mol), linewidth=1)

    for synth in synthetics : 
        AX2 = syntspec(path+synth)
        ax[1].plot(AX2['wavelen'], AX2['flux'], label = "Synth : "+ synthetics.get(synth), linewidth = 1)

    # ax[0].plot(AX['wavelen'], AX['flux'], label="Synth : "+ filename, linewidth=1)
    # ax[2].plot(AX2['wavelen'], AX2['flux'], label = "Synth : atom + mols", color = "royalblue", linewidth = 0.8)

    # ax[2].set_xlim(k-5,k+5)
    ax[1].set_xlim(k-10,k+10)
    ax[0].set_xlim(k-10,k+10)

    for ax in ax : 
        ax.set_xlabel("Longueur d'onde (Å)", fontsize  = 12)
        ax.set_ylabel("Flux normalisé", fontsize  = 12)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major')
        ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor')
        ax.set_ylim(0.2, 1.3)
        ax.axvline(x=k, ymin=0., ymax=0.9, color='black',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
        ax.text(k, 1.2, s = i, color='black', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre
        ax.legend(loc = "lower left", fontsize = 10)


def lines_analyse(lines, path, synthetics, molecules, stardata) :
    for i in list(lines.keys()):
        for k in lines.get(i):
            f = plt.figure(figsize=(15,8), dpi = 400)
            gs = f.add_gridspec(2, hspace=0.2)
            ax = gs.subplots(sharex=False, sharey=True)
            if k < 18500 : 
                zoom_lines_analyse(ax, path,synthetics, molecules, stardata, "", "h", k, i)

            elif k < 22500 : 
                zoom_lines_analyse(ax, path, synthetics,molecules, stardata, "2", "k", k, i)
            
            else :  
                zoom_lines_analyse(ax, path,synthetics, molecules, stardata, "3", "k", k, i)

            # plt.savefig("/Users/margauxvandererven/Documents/unif2023-2024/spectre_IR/output/"+stardata.get("starname")+"/"+i+"_"+str(k)+".png")


def zoom_lines_analyse_simple(ax, path,lines_mol, synthetics, stardata, m, k, i) : 
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)),stardata.get("flux_" + m),k)

    ax.plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='lightgray', label="Spectre observé : " + stardata.get("starname"))

    for synth in synthetics : 
        AX2 = syntspec(path+synth)
        ax.plot(AX2['wavelen'], AX2['flux'], label = "Synth : "+ synth, linewidth = 0.8)

    ax.set_xlim(k-10,k+10)
    
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major')
    ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor')
    ax.set_ylim(0.2, 1.3)
    ax.axvline(x=k, ymin=0., ymax=0.9, color='gray',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
    ax.text(k, 1.2, s = i, color='gray', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre
    ax.legend(loc = "lower left", fontsize = 6)

    for d in list(lines_mol.keys()):
        t = lines_mol.get(d)
        ax.axvline(x=t, ymin=0., ymax=0.9, color='gray',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
        # ax.text(t, 1.2, s = d, color='gray', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre


def lines_analyse_simple(lines, path, lines_mol, synthetics, stardata) :
    for i in list(lines.keys()):
        for k in lines.get(i):
            f = plt.figure(figsize=(20,8), dpi = 400)
            gs = f.add_gridspec(1, hspace=0.2)
            ax = gs.subplots(sharex=False, sharey=True)
            if k < 18500 : 
                zoom_lines_analyse_simple(ax, path, lines_mol,synthetics, stardata, "h", k, i)
            elif k < 22500 : 
                zoom_lines_analyse_simple(ax, path, lines_mol, synthetics, stardata, "k", k, i) 
            else :  
                zoom_lines_analyse_simple(ax, path, lines_mol,synthetics, stardata, "k", k, i)