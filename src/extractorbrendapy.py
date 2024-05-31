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

from brendapy import BrendaParser

DEFAULT_FILE_NAME = 'brenda_2023_1.txt'
# =============================================================================

# Extracts data from Brenda using Brendapy as parser.
# To create datasets for ML training

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
    return 'extractorbrendapy' + formatagedate + '.json'


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


def find_shared_key_p_ld(d_index : dict, d_kinetic : dict, p_cine : str,
                         v_key_p_list_dict) -> dict:
    """
    Dictionary that gathers the indexes of substrates that are present for the
    desired parameters (stored as a dict list in brenda)

    Ex: Gives the list of substrate indexes that are present for KM and TN
    TN if these are the two types of parameters requested by the user
    {'diacetyl' : {'KM': [1,2,3], 'TN' : [4,5,6]}}
    
    Index not taken into account when data value equals : more = ?

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
            if str(d_kinetic[i_subst][v_key_p_list_dict]) != 'more = ?':
                if not (str(d_kinetic[i_subst][v_key_p_list_dict]) in d_index):
                    d_index[str(d_kinetic[i_subst][v_key_p_list_dict])] = {p_cine : [i_subst]}
                elif not(p_cine in d_index[str(d_kinetic[i_subst][v_key_p_list_dict])]):
                    d_index[str(d_kinetic[i_subst][v_key_p_list_dict])].update({p_cine : [i_subst]})
                else:
                    d_index[str(d_kinetic[i_subst][v_key_p_list_dict])][p_cine].append(i_subst)
        except KeyError:
            pass
            logging.warning('Exception of key error')
    return d_index


def d_comment_each_kinetic(d_index : dict, d_i_substr : dict,
                           dict_proteins : dict) -> dict:
    """
    Associates all parameter comments with its index

    Parameters
    ----------
    d_index : dict
        
    d_i_substr : dict
        dictionary for each desired parameter with the list of indexes  in the
        brenda database for that substrate
    dict_proteins : OrderedDict
        Brenda is dictionary database

    Returns
    -------
    d_index : dict(dict)
        Dictionary whose key is the desired parameters, whose value is a
        dictionary with the index and its comment in the Brenda is database.
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
                # Problem with certain comments where the value is '-'.
                # In this case, the comment is ignored and not added to the database.
    return d_index


def find_keys_with_similar_values(main_dict: dict) -> list[dict]:
    """
    This function iterates through the nested structure of the input 
    and identifies keys that share common values.
    It returns a list of dictionaries, where each dictionary represents
    a group of keys with the same values.

    Finds keys in a nested dictionary with similar values

    Parameters
    ----------
    main_dict : dict
        Dictionary with imbricated sub-dictionaries

    Returns
    -------
    l_keys : List[Dict]
        Dictionary list containing indexes of parameters with similar comments

    Example
    -------
    Taking the output of the d_comment_each_kinetic function as an example
    Gives us :
    [{'TN': '16', 'KM': '20'}]

    """
    l_keys = []
    inverse_dict = defaultdict(list)

    for key, sub_dict in main_dict.items():
        for sub_key, comment in sub_dict.items():
            inverse_dict[comment].append((key, sub_key))

    # Recovers keys with common values
    keys_with_common_comment = {k: v for k, v in inverse_dict.items() if len(v) > 1}
    for v in keys_with_common_comment.values():
        l_keys.append(dict(v))
    return l_keys


def create_subdict_json(d_result, d_p_setting : dict, dict_proteins : dict,
                        i_sub_d_brenda, p_kinetic):
    """
    Put information extracted from Brenda in JSON format

    This function extracts data from BRENDAs stored in dict_proteins by
    selecting the parameter values specified in specified in d_p_setting and
    adds them to d_result.

    If a required key is missing in dict_proteins, it is ignored.

    Parameters
    ----------
    d_result : dict
        Results dictionary for data extracted from brenda
    d_p_setting : dict
        Dictionary of configuration parameters.
    dict_proteins : OrderedDict
        Brenda database in dictionary format
    i_sub_d_brenda :
        Brenda is sub-dictionary index
    p_kinetic :
        Parameter which is stored as list(dict) like KM or TN ...

    Returns
    -------
    d_result : dict
        Dict with data brenda
        Dictionary with data extracted from brenda

    """
    for parameter in d_p_setting['p_str']:
        d_result[str(parameter)] = dict_proteins[parameter]
    for parameter_k in d_p_setting['key_p_list_dict']:
        try:
            value_parameter = dict_proteins[p_kinetic][i_sub_d_brenda][parameter_k]
            # if value_parameter == 'more = ?':
            #     logging.warning('value : more = ? not accepted')
            # else:
            d_result[str(p_kinetic + '_' + parameter_k)] = value_parameter
        except KeyError:
            pass
            logging.warning('Exception of key error. \
                            Because, a comment with no value is ignored')
            # Problem with certain comments where the value is '-'.
            # In this case, the comment is ignored and not added to the database.
    return d_result


def check_parameter_values(d_p_setting : dict, dict_proteins):
    """
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


def k_subdict_parameter(cine):
    '''
    Depending on the parameter, which key to use to retrieve the value.

    Parameters
    ----------
    cine : str
        parameter.

    Returns
    -------
    str
        'substrate' or 'data'

    '''
    if cine in ['KM', 'KKM', 'KI', 'TN', 'IC50']:
        return 'substrate'
    return 'data'


def find_shared_substrate_index(para_list_dict : list, protein_data : dict) -> dict:
    """
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
    d_index_subkey : dict
        Dictionary for each desired parameter, which itself contains a 
        dictionary with substrates as key and a dictionary value for each 
        desired parameter with the list of indexes in brenda is database for
        that substrate.

    """
    d_index_subst = {}

    for cine in para_list_dict: #d_p_setting['p_list_dict']
        #substrat test
        # d_index_subkey = find_shared_substrate(d_index_subst,
        #                                       protein_data[cine], cine)
        # #no substrat -> data
        # if not d_index_subkey:
        #     d_index_subkey = find_shared_data(d_index_subst,protein_data[cine], cine)
        v_key_p_list_dict = k_subdict_parameter(cine)
        d_index_subkey = find_shared_key_p_ld(d_index_subst, protein_data[cine],
                                              cine,v_key_p_list_dict)
    return d_index_subkey


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
    results = []

    for ec_number in list_ec:
        for dict_proteins in BRENDA_PARSER.get_proteins(ec_number).values():
            if check_parameter_values(d_p_setting, dict_proteins.data):
                d_index_subst = find_shared_substrate_index(d_p_setting['p_list_dict'],
                                                            dict_proteins.data)

                for d_i_substr in d_index_subst.values():
                    d={}
                    d_index_comment = {}
                    for p_k, l_i_subst in d_i_substr.items():
                        #Those with a substrate
                        #All values are inserted directly into result
                        if len(l_i_subst) == 1:
                            d = create_subdict_json(d, d_p_setting,
                                                    dict_proteins.data,
                                                    l_i_subst[0], p_k)

                        #Those with several substrates for one parameter
                        #The locations of the values to be inserted in result
                        #for this parameter for this protein are temporarily
                        #stored in a dict / list so that they can be sorted,
                        #and then inserted separately.
                        #If we have several comments for the same substrate, 
                        #we'll put each comment in a different dictionary_json,
                        #but each time with the same general parameters,
                        #which don't change.
                        if len(l_i_subst) > 1:
                            d_index_comment = d_comment_each_kinetic(d_index_comment,
                                                                      d_i_substr,
                                                                      dict_proteins.data)
                            l_index_comment = find_keys_with_similar_values(d_index_comment)

                    #Those who have several comments, which can be different
                    #for the same parameter and the same substrate.
                    #But there may be comments that are the same for different
                    #parameters, in which case they are grouped together.
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

class DataSetBrenda:
    def __init__(self, list_paramaters : list, path_file_databrenda : str,
                 namefile : str, list_ec : list = []):

        if namefile:
            self.namefile = namefile
        else:
            self.namefile = name_new_file_created()

        self.d_parameter_setting = parameter_sorting(list_paramaters)
        self.path_set_brend = path_file_databrenda + self.namefile
        self.path_filebrenda = path_file_databrenda + DEFAULT_FILE_NAME

        self.list_ec = list_ec

    def get_cinetique_parameter(self):
        return self.d_parameter_setting
    def get_path_set_brend(self):
        return self.path_set_brend
    def get_namefile(self):
        return self.namefile

    def get_list_ec(self, Brendadata):
        if self.list_ec:
            return self.list_ec
        else:
            self.list_ec = Brendadata.keys()
            return self.list_ec

    def run(self):
        Brendadata = BrendaParser(self.path_filebrenda)
        create_file_json(self.get_path_set_brend(),
                         data_brenda(Brendadata, self.get_list_ec(Brendadata),
                                     self.get_cinetique_parameter()))