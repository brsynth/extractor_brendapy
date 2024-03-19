#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
import json
import logging

from collections import defaultdict
from datetime import datetime

# =============================================================================

# Creation une base de donnees qui servira comme jeux d'entrainement pour
# du machine learning

# =============================================================================

def name_new_file_created() -> str:
    """
    Gives the name of the file heure et date

    Returns
    -------
    str
        Names of the file to be created.

    """
    date_time = datetime.now()
    formatagedate = date_time.strftime('-%Y-%m-%d-%H-%M-%S')
    return 'setbrenda' + formatagedate + '.json'


def is_parameter_values(list_p : list, dict_proteins_data : dict) -> bool:
    """
    Checks that all the parameters in the list are in dict_proteins.data
    So if they have a value in this case, return TRUE

    Parameters
    ----------
    list_p : list
        list of parameters to be checked
    dict_proteins_data : OrderedDict
        Dictionary with all the parameters known to Brenda

    Returns
    -------
    bool
        return True or False.

    """
    for k_parameter in list_p:
        if not (k_parameter in dict_proteins_data):
            return False
        if not dict_proteins_data[str(k_parameter)]:
            return False
    return True


def find_shared_substrate(d_index : dict, d_kinetic : dict, p_cine : str) -> dict:
    """
    Dictionary that gathers the indexes of substrates that are present for the
    desired parameters (stored as a dict list in brenda)

    Ex: Gives the list of substrate indexes that are present for KM and TN
    TN if these are the two types of parameters requested by the user
    {'diacetyl' : {'KM': [1,2,3], 'TN' : [4,5,6]}}

    Parameters
    ----------
    d_index : dict
        dictionary with substrate as key and dictionary as value for each
        desired parameter with the list of indexes in brenda is database for
        that substrate.
    d_kinetic : dict
        brenda database for dictionary parameters
    p_cine : kinetic parameter

    Returns
    -------
    d_index : dict(dict)
        dictionary with substrate as key and dictionary as value for each
        desired parameter with the list of indexes in brenda is database for
        that substrate.

    """
    # index dictionary for different substrates
    for i_subst in range(len(d_kinetic)):
        try:
            if not (str(d_kinetic[i_subst]['substrate']) in d_index):
                d_index[str(d_kinetic[i_subst]['substrate'])] = {p_cine : [i_subst]}
            elif not(p_cine in d_index[str(d_kinetic[i_subst]['substrate'])]):
                d_index[str(d_kinetic[i_subst]['substrate'])].update({p_cine : [i_subst]})
            else:
                d_index[str(d_kinetic[i_subst]['substrate'])][p_cine].append(i_subst)
        except KeyError:
            pass
            logging.warning('Exception of key error')
    return d_index


def d_comment_each_kinetic(d_index : dict, d_i_substr : dict,
                           dict_proteins : dict) -> dict:
    """
    Associe tous les commentaires du parametre a son index

    Parameters
    ----------
    d_index : dict
        
    d_i_substr : dict
        dictionnaire pour chaque parametre souhaite avec la liste des index 
        dans la base de donne de brenda pour ce substrat.
    dict_proteins : OrderedDict
        Base de donnes de Brenda sous forme de dictionnaire

    Returns
    -------
    d_index : dict(dict)
        Dictionnaire Ayant pour clef les parametres souhaite, pour valeur un
        dictionnaire avec l'index et son commentaire dans la base de donne de
        Brenda.
        Ex: {
            'TN': {'16': 'pH 7.0, 25°C, mutant N107D <17>'},
             'KM': {'20': 'pH 7.0, 25°C, mutant N107D <17>',
                    '21': 'pH 7.0, 25°C, mutant N107L <17>'}}
    """
    for kinetic, l_index in d_i_substr.items():
        for index in l_index:
            try:
                if not (kinetic in d_index):
                    d_index[kinetic] = {str(index) : dict_proteins[kinetic][index]['comment']}
                elif not(index in d_index[kinetic]):
                    d_index[kinetic].update({str(index) : dict_proteins[kinetic][index]['comment']})
            except KeyError:
                pass
                logging.warning('Exception of key error. \
                                Because, a comment with no value is ignored')
                # Probleme avec certain commentaire ou la valeur est '-'
                # Dans ce cas, le commentaire est ignore donc pas ajoute a la 
                # base de donnees
    return d_index


def find_keys_with_similar_values(main_dict: dict) -> list[dict]:
    """
    Cette fonction itère à travers la structure imbriquée du dictionnaire 
    d'entrée et identifie les clés qui partagent des valeurs communes.
    Elle renvoie une liste de dictionnaires, où chaque dictionnaire représente
    un groupe de clés ayant les mêmes valeurs.

    Recherche les clés dans un dictionnaire imbrique ayant des valeurs
    similaires.

    This function iterates through the nested structure of the input 
    and identifies keys that share common values.
    It returns a list of dictionaries, where each dictionary represents
    a group of keys with the same values.

    Finds keys in a nested dictionary with similar values

    Parameters
    ----------
    main_dict : dict
        Dictionnaire avec des sous-dictionnaire imbriques.

    Returns
    -------
    l_keys : List[Dict]
        Liste de dictionnaire contenant les index des parametre qui ont des
        commentaires similaires.

    Exemples
    --------
    En reprenant l'exemple de sortie de la fonction d_comment_each_kinetic
    [{'TN': '16', 'KM': '20'}]

    """
    l_keys = []
    inverse_dict = defaultdict(list)

    for key, sub_dict in main_dict.items():
        for sub_key, comment in sub_dict.items():
            inverse_dict[comment].append((key, sub_key))

    # Recup les clés ayant des valeurs communes
    keys_with_common_comment = {k: v for k, v in inverse_dict.items() if len(v) > 1}
    for v in keys_with_common_comment.values():
        l_keys.append(dict(v))
    return l_keys


def create_subdict_json(d_result, d_p_setting : dict, dict_proteins : dict,
                        i_sub_d_brenda, p_kinetic):
    """
    Mets les informations extrait de Brenda au format JSON

    Cette fonction extrait les données des BRENDA qui sont stockes dans 
    dict_proteins en selectionnent les valeurs des parametres qui sont 
    specifies dans d_p_setting et les ajoute a d_result.

    Si une key necessaire est manquante dans dict_proteins, elle est ignoree.

    Put information extracted from Brenda in JSON format

    This function extracts data from BRENDAs stored in 
    dict_proteins by selecting the parameter values specified in 
    specified in d_p_setting and adds them to d_result.

    If a required key is missing in dict_proteins, it is ignored.

    Parameters
    ----------
    d_result : dict
        Dictionnaire des resultats des données extraites de brenda
    d_p_setting : dict
        Dictionary of configuration parameters.
    dict_proteins : OrderedDict
        Base de donnes de Brenda sous forme de dictionnaire.
    i_sub_d_brenda : TYPE
        index sous dictionnaire de Brenda.
    p_kinetic : TYPE
        parametre qui est stocke comme list(dict) comme KM ou TN ...

    Returns
    -------
    d_result : dict
        Dict with data brenda
        Dictionnaire avec les donnees extraites de brenda.

    """
    #mettre **kwarg
    for parameter in d_p_setting['p_str']:
        d_result[str(parameter)] = dict_proteins[parameter]
    for parameter_k in d_p_setting['key_p_list_dict']:
        try:
            value_parameter = dict_proteins[p_kinetic][i_sub_d_brenda][parameter_k]
            d_result[str(p_kinetic + '_' + parameter_k)] = value_parameter
        except KeyError:
            pass
            logging.warning('Exception of key error. \
                            Because, a comment with no value is ignored')
            # Probleme avec certain commentaire ou la valeur est '-'
            # Dans ce cas, le commentaire est ignore donc pas ajoute a la 
            # base de donnees
    return d_result


def check_parameter_values(d_p_setting : dict, dict_proteins):
    """
    Verifie la presence de valeurs dans la base de donnees pour les parametres
    souhaite.
    Checks the presence of values in the database for parameters wishes

    Parameters
    ----------
    d_p_setting : dict
        Dictionary of configuration parameters.
    dict_proteins : Orderdict
        dict_proteins.data.

    Returns
    -------
    Bool
        True or False.

    """
    get_params = [d_p_setting['p_str'], d_p_setting['p_list_dict']]
    return all(is_parameter_values(param, dict_proteins) for param in get_params)


def find_shared_substrate_index(para_list_dict : list, protein_data : dict) -> dict:
    """
    Lance la fonction find shared substrate pour tous les parametre souhaite
    qui sont dans la liste d_p_setting['p_list_dict'].

    Calls the find shared substrate function for all desired parameters
    in the d_p_setting['p_list_dict'] list.

    Parameters
    ----------
    para_list_dict : list
        list of parametres /d_p_setting['p_list_dict']
    protein_data : dict
        brenda protein data

    Returns
    -------
    d_index_subst : dict
        Dictionnaire pour chaque parametre souhaite, qui contient lui meme un 
        dictionnaire ayant en clef les substrat et en valeur un dictionnaire
        pour chaque parametre souhaite avec la liste des index dans la base de
        donne de brenda pour ce substrat.
        Dictionary for each desired parameter, which itself contains a 
        dictionary with substrates as key and a dictionary value for each 
        desired parameter with the list of indexes in brenda is database for
        that substrate.

    """
    d_index_subst = {}
    for cine in para_list_dict: #d_p_setting['p_list_dict']
        # print(cine, 'prot', protein_data)
        # print('prot data', protein_data[cine])
        d_index_subst = find_shared_substrate(d_index_subst,
                                              protein_data[cine], cine)
        # print('la 2', d_index_subst)
    return d_index_subst


def pre_subdict_from_couple(d_p_setting : dict, protein_data : dict,
                            couple : dict) -> dict:
    """
    creates subdictionary for parameters with multiple values linked by the 
    same comment

    Parameters
    ----------
    d_p_setting : dict
        Dictionary of configuration parameters.
    protein_data : TYPE
        brenda protein data.
    couple : dict
        comment with index.

    Returns
    -------
    dict
        dictionary in json format for identical comments on parameters with 
        multiple values.

    """
    d = {}
    for p_kine in d_p_setting['p_list_dict']:
        try:
            index_comment = int(couple[p_kine])
            d = create_subdict_json(d, d_p_setting, protein_data,
                                    index_comment, p_kine)
        except KeyError:
            pass
            logging.warning('Exception of key error.')
    return d


def data_brenda(BRENDA_PARSER, list_ec : list, d_p_setting : dict) -> list[dict]:
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

    #Pour cela nous allons cree une fonction qui trier les parametre en 3 classses:
        #Ceux qui peuvent etre recuperer direct
        #Ceux qui possede plusieurs valeurs -> list de dictionnaire
        #Ceux qui sont des parametre dans les dictionnaires : ex: substrate, comment, value
    """
    results = []

    for ec_number in list_ec:
        for dict_proteins in BRENDA_PARSER.get_proteins(ec_number).values():
            # print(dict_proteins.data)
            if check_parameter_values(d_p_setting, dict_proteins.data):
                d_index_subst = find_shared_substrate_index(d_p_setting['p_list_dict'],
                                                            dict_proteins.data)

                for d_i_substr in d_index_subst.values():
                    d={}
                    d_index_comment = {}
                    for p_k, l_i_subst in d_i_substr.items():
                        #Ceux qui possedent un substrat
                        #Tous les valeurs sont directement insere dans dans result
                        if len(l_i_subst) == 1:
                            d = create_subdict_json(d, d_p_setting,
                                                    dict_proteins.data,
                                                    l_i_subst[0], p_k)

                        #Ceux qui possedent plusieurs substrats pour un ce parametre
                        #Les localisations des valeurs a insere dans result pour
                        #ce parametre pour cette proteine sont temporairement stocke
                        #dans un dict / list pour pouvoir les triers, pour ensuite
                        #les inseres separements
                        #Siu nous avons plusieurs commentaire pour un meme substrat
                        #nous allons mettre chaque commentaire dans des dictionnaire_json
                        #different mais a chaque fois avec les meme parametre
                        #generaux qui eux ne change pas.
                        if len(l_i_subst) > 1:
                            d_index_comment = d_comment_each_kinetic(d_index_comment,
                                                                      d_i_substr,
                                                                      dict_proteins.data)
                            l_index_comment = find_keys_with_similar_values(d_index_comment)
                            # print(d_index_comment, l_index_comment)

                    #Ceux qui possedent plusieurs commentaire, qui peuvent 
                    #etre different pour le meme parametre et le meme substrat.
                    #Mais il peut y avoir des commentaires qui seront les meme
                    #pour des parametres different dans ce cas ils sont mis ensemble
                    if d_index_comment:
                        for couple in l_index_comment:
                            d = pre_subdict_from_couple(d_p_setting,
                                                        dict_proteins.data,
                                                        couple)
                            results.append(d)
                    else:
                        results.append(d)

    return results


def create_file_json(path_json : str, data : list[dict]):
    """
    Creates the json file and writes information on the various proteins with
    the parameters selected inside

    Parameters
    ----------
    path_json : strp_d_kinetic
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
    selection of elements common to both lists without redundancy

    Parameters
    ----------
    list1 : list
        list of parametres.
    list2 : list
        list of parametres.

    Returns
    -------
    list
        list containing only parameters common to both lists

    """
    return list(set(list1) & set(list2))


def parameter_sorting(list_parameter : list) -> dict:
    """
    Sorts parameters by type in Brenda's database

    Parameters
    ----------
    list_parameter : list
        list of user-defined parameters

    Returns
    -------
    Dict
        Classification of user-defined parameters

    """
    all_parameter = {'p_str' : ["ec", "uniprot", "organism", "ID"],
                 'p_list_dict' : ['KM', 'KKM', 'KI', 'TN', 'IC50', "ref", "TS",
                                  "SY", "SU", "ST", "SP", "SA", "PU", "NSP",
                                  "MW", "LO", "GI", "IN", "CL", "CF", "AP"],
                 'key_p_list_dict' : ['substrate', 'value', 'comment', 'units',
                                      'refs', 'data', 'chebi'],
                 'p_set' : ["tissues", "SN", "RT", "RN", "RE"]}
    d_parameter_setting = {'p_str' : commun_lists(list_parameter,
                                                  all_parameter['p_str']),
                 'p_list_dict' : commun_lists(list_parameter,
                                              all_parameter['p_list_dict']),
                 'key_p_list_dict' : commun_lists(list_parameter,
                                                  all_parameter['key_p_list_dict']),
                 'p_set' : commun_lists(list_parameter,
                                        all_parameter['p_set'])}
    return d_parameter_setting
