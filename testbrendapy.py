#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from typing import Dict, List
import json

# from collections import defaultdict

# path_data_brenda = '/home/nparis/brenda_enzyme/'
# file_name_json = 'brenda_2023_1.json'
# file_name_txt = 'brenda_2023_1.txt'

from brendapy import BrendaParser, BrendaProtein
# from brendapy.console import console
# from brendapy.log import get_logger
# from brendapy.taxonomy import Taxonomy

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
# reorganisation du code et re-ecriture de la docstring qui n'est plus a jours!
# =============================================================================

def file_path_request(path_brenda : str) -> str:
    """
    fourni le chemin jsuqu'au fichier qui contient tous les data de brenda au
    format txt

    Parameters
    ----------
    path_brenda : str
        chemin jusqu'au dossier ou est situer le fichier txt.

    Returns
    -------
    str
        chemin avec le noms du fichier txt inclus.

    """
    file_name_txt = 'brenda_2023_1.txt'
    return str(path_brenda+file_name_txt)


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


def is_parameter_values(list_p : list, dict_proteins_data : dict) -> bool:
    """
    verifie que tous les parametre de la liste sont dans dict_proteins.data
    S'ils sont dedans retourne TRUE

    Parameters
    ----------
    list_p : list
        liste des parametre a verifier.
    dict_proteins_data : OrderedDict
        ictionnaire ou il y a tous les parametre connu dans Brenda.

    Returns
    -------
    bool
        reourne True ou False.

    """
    #TODO : concateniser
    for k_parameter in list_p:
        # n'accepte pas le if not ... and ... :
        if not (k_parameter in dict_proteins_data):
            return False
        if not dict_proteins_data[str(k_parameter)]:
            return False
    return True


def find_shared_substrate(d_index : dict, d_kinetic : dict) -> dict:
    """
    dictionnaire contenant les index des substrats qui sont partager pour les
    differrents parametre de cinetique demande

    Parameters
    ----------
    d_index : dict
        DESCRIPTION.
    d_kinetic : dict
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    # dictionnaire index pour les differents substrats
    for i_subst in range(len(d_kinetic)):
        try:
            if not (d_kinetic[i_subst]['substrate'] in d_index):
                d_index[d_kinetic[i_subst]['substrate']] = [i_subst]
            else:
                d_index[str(d_kinetic[i_subst]['substrate'])].append(i_subst)
        except KeyError:
            pass
    return d_index


def data_brenda(list_ec : list, d_p_setting : dict) -> List[Dict]:
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
    d_p_setting : dict
        

    Returns
    -------
    List[Dict]
        Parameter dictionary list for each proteins.

    """
    """
    #Dans un premier temps on va mettre dans le dictionnaire les parametres
    #ou nous y avons acces directement
    
    #Pour cela nous allons cree une fonction qui trier les parametre en 3 classses:
        #Ceux qui peuvent etre recuperer direct
        #Ceux qui possede plusieurs valeurs -> list de dictionnaire
        #Ceux qui sont des parametre dans les dictionnaires : ex: substrate, comment, value
    """

    results = []
    for ec_number in list_ec:
        for dict_proteins in BRENDA_PARSER.get_proteins(ec_number).values():
            #verifie que tous les parametre de la liste sont dans dict_proteins.data
            #si ils sont dedans -> possede une valeur
            if is_parameter_values(d_p_setting['p_primaire'], dict_proteins.data) and is_parameter_values(d_p_setting['p_kinetic'], dict_proteins.data):
                # d = {}
                d_tempo = {}
                d_index_subst = {}
                for cine in d_p_setting['p_kinetic']:
                    #mettre toutes les parametre de cinetique ensemble pour les reorganise par substrat
                    #s'ils ont le meme substrat on mets les information km et tn ensemble sinon on les mets dans des dictionnaire differents
                    #pour eviter d'avoir des doublons
                    d_tempo[str(cine)] = dict_proteins.data[cine]
                    
                    d_index_subst = find_shared_substrate(d_index_subst, dict_proteins.data[cine])
                    # print(d_index_subst)
                #Dans un pemiere temps nousallons nous occupe des subst qui
                #ne sont present qu'une fois
                # TODO : les substrats presents plusieurs fois, ayant des comment differents
                for substr in d_index_subst.keys():
                    if len(d_index_subst[substr]) == len(d_p_setting['p_kinetic']):
                        d={}
                        for i in range(len(d_p_setting['p_kinetic'])):
                            # str du parametre de cine : d_p_setting['p_kinetic'][i])
                            #index dans le dictionnaire de cinetique : d_index_subst[substr][i])
                            # print(d_tempo[str(d_p_setting['p_kinetic'][i])][d_index_subst[substr][i]])
                            d[str(d_p_setting['p_kinetic'][i])] = d_tempo[str(d_p_setting['p_kinetic'][i])][d_index_subst[substr][i]]
                        for p in d_p_setting['p_primaire']:
                            d[str(p)] = dict_proteins.data[p]
                        results.append(d)
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


def commun_lists(list1 : list, list2 : list) -> list:
    """
    selection des elements commun entre les deux listes sans redondance

    Parameters
    ----------
    list1 : list
        liste de parametres.
    list2 : list
        liste de parametres.

    Returns
    -------
    list
        liste contenant uniquement les parametres communsentre les deux listes.

    """
    return list(set(list1) & set(list2))


#Verifier que les parameter mis en argument appartient bien au dict
def parameter_sorting(list_parameter : list) -> Dict:
    """
    Trie les paramtres en fonction de leur placement dans Brenda pour les recuperer

    Parameters
    ----------
    list_parameter : list
        DESCRIPTION.

    Returns
    -------
    Dict
        DESCRIPTION.

    """
    all_parameter = {'p_primaire' : ["ec", "uniprot", "organism"],
                 'p_kinetic' : ['KM', 'KKM', 'KI', 'TN', 'IC50'],
                 'p_d_kinetic' : ['substrate', 'value', 'comment', 'units', 'refs']}
    d_parameter_setting = {'p_primaire' : commun_lists(list_parameter, all_parameter['p_primaire']),
                 'p_kinetic' : commun_lists(list_parameter, all_parameter['p_kinetic']),
                 'p_d_kinetic' : commun_lists(list_parameter, all_parameter['p_d_kinetic'])}
    return d_parameter_setting


class DataSetBrenda:
    def __init__(self, list_paramater : list, path_data_brenda : str):
        self.d_parameter_setting = parameter_sorting(list_paramater)
        #noms du fichier avec la date et heure de creation ...
        self.path_set_brend = path_data_brenda + name_new_file_created(self.d_parameter_setting['p_kinetic'][0])

    def get_cinetique_parameter(self):
        return self.d_parameter_setting
    def get_path_set_brend(self):
        return self.path_set_brend

    #@property -> enlever () apres le run dans l'appel de fonction
    def run(self):
        create_file_json(self.get_path_set_brend(), data_brenda(list_all_ec_in_data(),
                                                              self.get_cinetique_parameter()))


# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

BRENDA_PARSER = BrendaParser(file_path_request('/home/nparis/brenda_enzyme/'))


# create_file_json(str(path_data_brenda+name_new_file_created('KM')),
#                  data_brenda(list_all_ec_in_data(),'KM'))

list_p = ["ec", "uniprot", "organism", "substrate", 'comment', 'KM', 'TN']
# DataSetBrenda(list_p, '/home/nparis/brenda_enzyme/').run()

d_parameter_setting = parameter_sorting(list_p)
brendaset = data_brenda(['1.1.1.10'], d_parameter_setting)
# create_file_json(str(path_data_brenda+name_new_file_created('KM')), data_brenda(['1.1.1.1'], 'KM'))

# print(BRENDA_PARSER.BRENDA_KEYS)


# Affichage Brendapy
# console.rule()
# console.print(brendaset)
# console.rule()