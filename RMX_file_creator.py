#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:55:58 2024

@author: nparis
"""
import json
import re

# Lire le fichier
# Lire la ligne data
# Ajouter au fichier susbtrat et produit


# path = '/home/nparis/brenda_enzyme/'
# input_file = 'test1.json'
# file_out = 'out1.json'

def molecule_sep(elements: str):
    """
    

    Parameters
    ----------
    elements : str
        DESCRIPTION.

    Returns
    -------
    result : TYPE
        DESCRIPTION.

    """
    eq = elements.split('+')
    result = []
    
    for elet in eq:
        elet = elet.strip()
        
        match = re.match(r'(\d*)\s*(.+)', elet)
        if match:
            coef = match.group(1)
            molecule = match.group(2)
            if coef == '':
                coef = 1
            else:
                coef = int(coef)

            result.append((coef, molecule))

    return result

# print(molecule_sep('monohexosylceramide + 2 ferrocytochrome b5 + O2 + 2 H+'))
# print(molecule_sep('ATP + H2O'))


def modif_file(path : str, input_file : str, file_out : str):
    """
    

    Parameters
    ----------
    path : str
        DESCRIPTION.
    input_file : str
        DESCRIPTION.
    file_out : str
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with open(path + input_file, "r") as file:
        data = json.load(file)
    for element in data:
        reaction_SP = element['SP_data']
        i_symbol_egale = reaction_SP.find("=")
        substrates = molecule_sep(reaction_SP[:i_symbol_egale-1])
        produits = molecule_sep(reaction_SP[i_symbol_egale+2:])
    
        element['substrates'] = substrates
        element['products'] = produits
    
    with open(path + file_out, "w", encoding = 'utf8') as file:
        json.dump(data, file, indent = 2, ensure_ascii=False)


class RMXData:
    def __init__(self, path, input_file, file_out):
        self.path = path
        self.input_file = input_file
        if file_out:
            self.file_out = file_out
        # else:
        #     self.file_out = faire une fonction qui rajoute RMX dans le noms du fichier d entre

    def get_path(self):
        return self.path
    def get_input_file(self):
        return self.input_file
    def get_file_out(self):
        return self.file_out

    def run(self):
        modif_file(self.get_path(), self.get_input_file(), self.get_file_out())