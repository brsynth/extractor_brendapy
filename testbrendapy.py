#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from typing import Dict, List
import json

path_data_brenda = '/home/nparis/brenda_enzyme/'
# file_name_json = 'brenda_2023_1.json'
file_name_txt = 'brenda_2023_1.txt'
name_new_file_json = 'setbranda.json'

from brendapy import BrendaParser, BrendaProtein
# from brendapy.console import console
# from brendapy.log import get_logger
# from brendapy.taxonomy import Taxonomy

BRENDA_PARSER = BrendaParser(str(path_data_brenda+file_name_txt))

# =============================================================================
# Parameters we want to retrieve from Brenda:
# EC number
# Km values
# Substrat
# Organisms
# Uniprot
# Commentary

# Enregistrer les donnees dans un fichier au format json
# =============================================================================

#TODO : liste des parametres
#Dict des parametre souhaite
# dict_para = {"EC" : None,'protein_id':None, 'organism':None, 'uniprot':None, 'KM':None}
# def create_dict_para() -> Dict:
#     pass

def list_all_ec_in_data() -> List:
    """
    Gives all known EC in brenda

    Returns
    -------
with open("/chemin/vers/le_fichier.json", "w") as f:
    json.dump(data, f)
    List
        list all EC number in Brenda.

    """
    return BRENDA_PARSER.keys()


def data_brenda(list_ec : list) -> List[Dict]:
    """
    List containing a dictionary for each protein with the parameters
    selected
    By selecting proteins with only the desired parameters
    Liste contenant un dictionnaire pour chaque proteine avec les parametres
    selectionnes
    En selectionnant les proteines qui ont uniquement les parametres souhaite

    BRENDA_PARSER.get_proteins returns Dict[int, BrendaProtein]

    Parameters
    ----------
    list_ec : list
        All EC number in Brenda.

    Returns
    -------
    List[Dict]
        Parameter dictionary list for each proteins.

    """
    results = []
    for ec_number in list_ec:
        for k, dict_proteins in BRENDA_PARSER.get_proteins(ec_number).items():
            #TODO: facto la condition et les parametres pour pouvoir les modif
            if dict_proteins.uniprot and dict_proteins.KM:
                for dict_KM in dict_proteins.KM:
                    try:
                        results.append({"EC" : ec_number, "protein_id": k,
                                        "uniprot": dict_proteins.uniprot,
                                        "organism": dict_proteins.organism,
                                        "substrate": dict_KM['substrate'],
                                        'KM' : dict_KM['value'],
                                        'comment': dict_KM['comment']})
                    except KeyError:
                        '''Parfois nous n'avons pas l'information pour la valeur
                        du Km et le noms du substrat'''
                        pass
    return results


# r = data_brenda(['1.1.1.1'])

# brendaset = data_brenda(list_all_ec_in_data())

# Affichage Brendapy
# console.rule()
# console.print(brendaset)
# console.rule()

# =============================================================================
# Ecrire les donnees dans un fichier JSON
# =============================================================================

def create_file_json(path_json : str, data : list[Dict]):
    """
    

    Parameters
    ----------
    path_json : str
        DESCRIPTION.
    data : list[Dict]
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with open(path_json, "w") as file:
        json.dump(data, file, indent=2)

create_file_json(str(path_data_brenda+name_new_file_json), data_brenda(list_all_ec_in_data()))