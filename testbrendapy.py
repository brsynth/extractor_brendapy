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
# from brendapy.console import console
# from brendapy.log import get_logger
# from brendapy.taxonomy import Taxonomy

BRENDA_PARSER = BrendaParser(str(path_data_brenda+file_name_txt))

# =============================================================================
# Parameters we want to retrieve from Brenda:
# EC number
# Km /Kcat values
# Substrat
# Organisms
# Uniprot
# Commentary

# Save data in a json format file

# List parametre souhaite
# reorganisation du code
# contacte la personne de brendapy
# =============================================================================

def name_new_file_created(cinetique_parameter : str) -> str:
    """
    Gives the name of the file according to the selected parameters in json
    format

    Parameters
    ----------
    cinetique_parameter : str
        type of kinetics parameter we're going to retrieve and store in the
        file.

    Returns
    -------
    str
        Names of the file to be created.

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

            a = getattr(dict_proteins, cine_parameter)

            if dict_proteins.uniprot and a:
                #dict_cine = soit dict_KM soit dict_TN = Kcat
                for dict_cine in a:
                    try:
                        results.append({"EC" : ec_number,
                                        "uniprot": dict_proteins.uniprot,
                                        "organism": dict_proteins.organism,
                                        "substrate": dict_cine['substrate'],
                                        str(cine_parameter) : dict_cine['value'],
                                        'comment': dict_cine['comment']})
                    except KeyError:
                        '''Sometimes we don't have information on the Km value
                        and substrate names'''
                        pass
    return results


def create_file_json(path_json : str, data : list[Dict]):
    """
    Creates the json file and writes information on the various proteins with
    the parameters selected inside

    Parameters
    ----------
    path_json : str
        path of the json file to create.
    data : list[Dict]
        Parameter dictionary list for each proteins.

    Returns
    -------
    Creation of a file with the requested names at the requested location in
    json format.

    """
    with open(path_json, "w", encoding = 'utf8') as file:
        json.dump(data, file, indent = 2, ensure_ascii=False)


class DataSetBrenda:
    def __init__(self, cinetique_parameter : str, path_data_brenda : str):
        #cinetique parameter = Km ou Kcat (=TN)
        self.cinetique_parameter = cinetique_parameter
        self.path_set_brend = path_data_brenda + name_new_file_created(cinetique_parameter)

    def get_cinetique_parameter(self):
        return self.cinetique_parameter
    def get_path_set_brend(self):
        return self.path_set_brend

    #@proprerty -> enlever () apres le run dans l'appel de fonction
    def run(self):
        create_file_json(self.get_path_set_brend(), data_brenda(list_all_ec_in_data(),
                                                             self.get_cinetique_parameter()))



# if __name__ == '__main__':

# create_file_json(str(path_data_brenda+name_new_file_created('KM')),
#                  data_brenda(list_all_ec_in_data(),'KM'))

DataSetBrenda('TN', path_data_brenda).run()

# brendaset = data_brenda(['1.1.1.1'], 'KM')
# r = data_brenda(['1.1.1.1'], 'TN')

# brendaset = data_brenda(list_all_ec_in_data())

# Affichage Brendapy
# console.rule()
# console.print(brendaset)
# console.rule()