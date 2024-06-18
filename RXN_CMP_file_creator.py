#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:55:58 2024

@author: nparis
"""
import json
import re
import requests
from rdkit import Chem
from zeep import Client
# Lire le fichier
# Lire la ligne data
# Ajouter au fichier susbtrat et produit


path = '/home/nparis/brenda_enzyme/'
input_file = 'test1.json'
file_out = 'out1.json'

# =============================================================================
# PARTIE RXN
# =============================================================================

def molecule_sep(elements: str):
    """
    Separates a string of chemical elements into a list of tuples containing 
    coefficients and molecules.
    
    Parameters
    ----------
    elements : str
        A string representing chemical elements and their coefficients
    
    Returns
    -------
    list of tuple
        A list of tuples, where each tuple contains an integer coefficient 
        and a string molecule. JSON tuples = list
    """
    eq = elements.split('+')
    result = []
    l_molecule = {}
    
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

            l_molecule[molecule] = ''
            result.append((coef, molecule))

    return result, l_molecule

# print(molecule_sep('monohexosylceramide + 2 ferrocytochrome b5 + O2 + 2 H+'))
# print(molecule_sep('ATP + H2O'))


def modif_file(path : str, input_file : str, file_out : str):
    """
    Modifies a JSON file by parsing chemical reactions and adding substrate 
    and product information.
    
    Parameters
    ----------
    path : str
        The directory path where the input and output files are located
    input_file : str
        The name of the input JSON file containing reaction data
    file_out : str
        The name of the output JSON file where modified data will be saved
    
    Returns
    -------
    This function does not return any value. It writes the modified data to
    the specified output file.
    """
    with open(path + input_file, "r") as file:
        data = json.load(file)
    CMP_data = {}
    for element in data:
        reaction_SP = element['SP_data']
        i_symbol_egale = reaction_SP.find("=")
        substrates = molecule_sep(reaction_SP[:i_symbol_egale-1])
        produits = molecule_sep(reaction_SP[i_symbol_egale+2:])
        CMP_data.update(substrates[1])
        CMP_data.update(produits[1])
        # print()
        # print(CMP_data)
        # re-extrait les ID a partir de sub et prd
        #ou
        # mol sep retourne result et la list des ID pour faire le fichier CMP
        #enregistre avec une chaine vide et c'est plus tard que je mets le smile
        element['substrates'] = substrates[0]
        element['products'] = produits[0]
    
    with open(path + file_out, "w", encoding = 'utf8') as file:
        json.dump(data, file, indent = 2, ensure_ascii=False)
    file_out2 = 'out2.json'
    with open(path + file_out2, "w", encoding = 'utf8') as file:
        json.dump(CMP_data, file, indent = 2, ensure_ascii=False)


def new_filename(file: str, name : str) -> str:
    """
    Adds 'RXN' to the filename.
    
    Parameters
    ----------
    file : str
        The name of the file to which 'RXN' will be added.
    
    Returns
    -------
    new_filename : str
        The new filename with 'RXN' added at the beginning.
    """
    new_file = f"{name}_{file}"
    return new_file

# =============================================================================
# PARTIE CMP
# =============================================================================

def get_pubchem_cid_from_name(protein):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{protein}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID'][0]
    return None

# p = "ATP"
# pubchem_cid = get_pubchem_cid_from_name(p)


def get_molfile_from_pubchem(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

# molfile_data = get_molfile_from_pubchem(pubchem_cid)


def molfile_to_smiles2(molfile):
    try:
        # print(molfile)
        mol = Chem.MolFromMolBlock(molfile)
        if mol is None:
            raise ValueError("nope premier etape de conversion en smile")
        smiles = Chem.MolToSmiles(mol)
        return smiles

    except Exception as e:
        print('le prbl est :', e)

# smiles = molfile_to_smiles2(molfile_data)
# print(smiles)


# =============================================================================
# PARTIE FUSION RXN / CMP
# =============================================================================
class RXNData:
    def __init__(self, path, input_file, file_out):
        self.path = path
        self.input_file = input_file
        if file_out:
            self.file_out = file_out
        else:
            self.file_out = new_filename(input_file, 'RXN')

    def get_path(self):
        return self.path
    def get_input_file(self):
        return self.input_file
    def get_file_out(self):
        return self.file_out

    def run(self):
        modif_file(self.get_path(), self.get_input_file(), self.get_file_out())

RXNData(path, input_file, file_out).run()