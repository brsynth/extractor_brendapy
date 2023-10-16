#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from typing import Dict, List
import json

# from collections import defaultdict

path_data_brenda = '/home/nparis/brenda_enzyme/'
# file_name_json = 'brenda_2023_1.json'
file_name_txt = 'brenda_2023_1.txt'

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


def name_new_file_created(cinetique_parameter : str) -> str:
    """
    

    Parameters
    ----------
    cinetique_parameter : str
        type de parametre de cinetique que nous allons recuperer et stocker
        dans le ficier.

    Returns
    -------
    str
        Noms du fichier a cree.

    """
    return 'setbrenda_' + str(cinetique_parameter) + '.json'


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
                #dict_cine = soit dict_KM soit dict_KKM = Kcat
                for dict_cine in a:
                    try:
                        results.append({"EC" : ec_number,
                                        "uniprot": dict_proteins.uniprot,
                                        "organism": dict_proteins.organism,
                                        "substrate": dict_cine['substrate'],
                                        str(cine_parameter) : dict_cine['value'],
                                        'comment': dict_cine['comment']})
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
    Creation d'un fichier au noms demande a l'emplacement demande au format json.

    """
    with open(path_json, "w") as file:
        json.dump(data, file, indent=2)

# create_file_json(str(path_data_brenda+name_new_file_created('KM')), data_brenda(list_all_ec_in_data(),'KM'))


class DataSetBrenda:
    def __init__(self, cinetique_parameter : str, path_data_brenda : str):
        self.cinetique_parameter = cinetique_parameter
        self.path_set_brend = path_data_brenda + name_new_file_created(cinetique_parameter)

    def get_cinetique_parameter(self):
        return self.cinetique_parameter
    def get_path_set_brend(self):
        return self.path_set_brend

    def run(self):
        create_file_json(self.get_path_set_brend(), data_brenda(list_all_ec_in_data(),
                                                             self.get_cinetique_parameter()))


# print(DataSetBrenda('KKM', path_data_brenda).get_path_set_brend())
DataSetBrenda('KKM', path_data_brenda).run()
