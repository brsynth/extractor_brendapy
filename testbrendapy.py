#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from typing import Dict, List
import json

from collections import defaultdict

path_data_brenda = '/home/nparis/brenda_enzyme/'
# file_name_json = 'brenda_2023_1.json'
file_name_txt = 'brenda_2023_1.txt'
name_new_file_json = 'setbranda.json'

from brendapy import BrendaParser, BrendaProtein
from brendapy.console import console
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

#Dict des parametre souhaite qui seront un argument de la classe -> kwarg


# def desired_parameter(ec_number : str, dict_proteins : dict, cine_parameter) -> Dict:
#     """
    

#     Returns
#     -------
#     Dict
#         DESCRIPTION.

#     """
#     if cine_parameter == 'KM':
#         a = dict_proteins.KM
#     else:
#         a = dict_proteins.KKM

#     if dict_proteins.uniprot and a:
#         for dict_cinetique in a:
#             try:
#                 dict_result = {"EC" : ec_number,
#                                 "uniprot": dict_proteins.uniprot,
#                                 "organism": dict_proteins.organism,
#                                 "substrate": dict_cinetique['substrate'],
#                                 str(cine_parameter) : dict_cinetique['value'],
#                                 'comment': dict_cinetique['comment']}
#             except KeyError:
#                 '''Parfois nous n'avons pas l'information pour la valeur
#                 du Km ou le KKM (=Kcat) et le noms du substrat'''
#                 pass
#     return dict_result


def list_all_ec_in_data() -> List:
    """
    Gives all known EC in brenda

    Returns
    -------
    List
        list all EC number in Brenda.

    """
    return BRENDA_PARSER.keys()


def data_brenda(list_ec : list, cine_parameter : str) -> List[Dict]:
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
    parameter : str

    Returns
    -------
    List[Dict]
        Parameter dictionary list for each proteins.

    """
    results = []
    for ec_number in list_ec:
        for k, dict_proteins in BRENDA_PARSER.get_proteins(ec_number).items():

            if cine_parameter == 'KM':
                a = dict_proteins.KM
            else:
                a = dict_proteins.KKM

            if dict_proteins.uniprot and a:
                for dict_KM in a:
                    try:
                        results.append({"EC" : ec_number,
                                        "uniprot": dict_proteins.uniprot,
                                        "organism": dict_proteins.organism,
                                        "substrate": dict_KM['substrate'],
                                        str(cine_parameter) : dict_KM['value'],
                                        'comment': dict_KM['comment']})
                    except KeyError:
                        '''Parfois nous n'avons pas l'information pour la valeur
                        du Km et le noms du substrat'''
                        pass
    return results


# brendaset = data_brenda(['1.1.1.1'], 'KM')
# r = data_brenda(['1.1.1.1'], 'KKM')

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

# create_file_json(str(path_data_brenda+name_new_file_json), data_brenda(list_all_ec_in_data()))



class DataSetBrenda():
    pass





