#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from collections import defaultdict
from typing import Dict, List

path_data_brenda = '/home/nparis/brenda_enzyme/'
file_name_json = 'brenda_2023_1.json'
file_name_txt = 'brenda_2023_1.txt'

from brendapy import BrendaParser, BrendaProtein
from brendapy.console import console
from brendapy.log import get_logger
from brendapy.taxonomy import Taxonomy

BRENDA_PARSER = BrendaParser(str(path_data_brenda+file_name_txt))

# =============================================================================
# Parameters we want to retrieve from Brenda:
# EC number
# Km values
# Substrat
# Organisms
# Uniprot
# Commentary
# =============================================================================

#TODO : liste des parametres
#Dict des parametre souhaite
# dict_set = {"EC" : None,'protein_id':None, 'organism':None, 'uniprot':None, 'KM':None}
# def create_dict_set() -> Dict:
#     pass

def list_all_EC_in_data() -> List:
    """
    Gives all known EC in brenda

    Returns
    -------
    List
        list all EC number in Brenda.

    """
    return BRENDA_PARSER.keys()


def data_brenda(list_EC : list) -> List[Dict]:
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
    list_EC : list
        All EC number in Brenda.

    Returns
    -------
    List[Dict]
        Parameter dictionary list for each proteins.

    """
    results = list()
    for EC_number in list_EC:
        for k, v in BRENDA_PARSER.get_proteins(EC_number).items():
            #TODO: facto la condition et les parametres pour pouvoir les modif
            if v.uniprot and v.KM:
                results.append({"EC" : EC_number, "protein_id": k,
                                "uniprot": v.uniprot,"organism": v.organism,
                                "KM": v.KM})
    return results


r = data_brenda(list_all_EC_in_data())

# Affichage Brendapy
# console.rule()
# console.print(results)
# console.rule()









































