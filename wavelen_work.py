import numpy as np

def syntspec(pathtofile):
    """
    Lecture des fichiers synthétiques
    """
    file = open(pathtofile)
    wavelen = []
    flux = []
    for line in file:
        i = line.split() 
        wavelen.append(float(i[0]))
        flux.append(float(i[1]))
    return {'wavelen' : wavelen, 'flux' : flux }


def redshift_wavelen(wavelen, v):
    """
    Fonction qui décalent les longueurs d'onde en fonction de redshift
    du à la vitesse radiale de l'étoile
    """
    redshift = []
    c = 299792458 #m/s
    for i in wavelen:
        z = v/c * i
        k = i - z
        redshift.append(k)
    return redshift   


def mediane_5(liste):
    liste_triee = sorted(liste, reverse=True)
    nb_elements_a_extraire = int(len(liste) * 0.05)
    plus_grands = liste_triee[:nb_elements_a_extraire]
    liste_triee2 = sorted(plus_grands)
    longueur = len(liste_triee2)
    if longueur % 2 != 0:
        mediane_index = longueur // 2
        mediane = liste_triee2[mediane_index]
    else:
        mediane_index_1 = longueur // 2 - 1
        mediane_index_2 = longueur // 2
        mediane = (liste_triee2[mediane_index_1] + liste_triee2[mediane_index_2]) / 2
    
    return mediane



def zoom_syntspec(path, filename, wavelength, taille):
    """
    Zoom sur le spectre synthétique pour une longueur d'onde donnée
    """
    wavelenB = syntspec(path+filename)['wavelen']
    fluxB = syntspec(path+filename)['flux']
    h = []
    for j in wavelenB:
        if 0 < j-wavelength:
            h.append(j)
    centre_synt = wavelenB.index(min(h))
    length = int(taille/0.03)
    
    return {"synt_wavelen" : wavelenB[centre_synt-length:centre_synt+length], "synt_flux" : fluxB[centre_synt-length:centre_synt+length]}


def zoom_syntspec2(path, filename, z_wavelen, flux, wavelength):
    """
    Zoom sur le spectre synthétique pour une longueur d'onde donnée, superposé au spectre observé
    """
    wavelenB = syntspec(path+filename)['wavelen']
    fluxB = syntspec(path+filename)['flux']
    k = []
    h = []
    for i in z_wavelen:
        if 0 < i-wavelength:
            k.append(i)
    centre = z_wavelen.index(min(k))

    for j in wavelenB:
        if 0 < j-wavelength:
            h.append(j)
    centre_synt = wavelenB.index(min(h))

    mediane = mediane_5(flux[centre-160:centre+160])

    flux_normalised = []
    for m in flux[centre-160:centre+160] :
        flux_normalised.append(m/mediane)

    return {"z_wavelen" : z_wavelen[centre-160:centre+160], "flux_normalised" : flux_normalised, 
            "synt_wavelen" : wavelenB[centre_synt-2000:centre_synt+2000], "synt_flux" : fluxB[centre_synt-2000:centre_synt+2000]}


def get_nearest(x_query, x_vals):
    return np.abs(x_vals - x_query).argmin()


def normalisation(z_wavelen, flux, wavelength):
    """
    Normalisation du flux observé sur interval de 20 Å. 
    Le flux est divisé par la médiane des 5% plus grandes valeurs sur cet interval.
    """
    # k = []
    h = []
    # for i in z_wavelen:
    #     if 0 < i-wavelength:
    #         k.append(i)
    # centre = z_wavelen.index(min(k))

    centre = get_nearest(wavelength, np.array(z_wavelen))

    mediane = mediane_5(flux[centre-160:centre+160])

    flux_normalised = []
    for h in flux[centre-160:centre+160] :
        flux_normalised.append(h/mediane)
    # flux_normalised = np.array()

    tolerance = 0.05
    index_closest_with_tolerance = next((i for i, x in enumerate(flux_normalised[160:]) if abs(x - 1) <= tolerance), None)

    # taille = z_wavelen[centre+index_closest_with_tolerance] - z_wavelen[centre]
    
    return {"z_wavelen" : z_wavelen[centre-160:centre+160], "flux_normalised" : flux_normalised, 
            # "taille" : taille,
             "index_closest_with_tolerance" : index_closest_with_tolerance}