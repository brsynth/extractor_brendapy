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
import hashlib
import logging
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
        A list of tuples, whimport loggingere each tuple contains an integer coefficient 
        and a string molecule. JSON tuples = list
    """
    #Pre traitement
    elements = re.sub(r'\|.*?\|', '', elements).strip()
    # #Laisser l'espace sinon ne prend pas le + du H+ par exemple
    eq = elements.split(' + ')
    molecules = {}
    result = []
    rever = '?'

    for elet in eq:
        elet = elet.strip()
        
        if elet == '?':
            return None

        if '{r}' in elet:
            rever = 'r'
            elet = elet.replace('{r}', '').strip()
        elif '{ir}' in elet:
            rever = 'ir'
            elet = elet.replace('{ir}', '').strip()

        match = re.match(r'(\d+)\s+(.+)', elet)
        if match:
            coef = int(match.group(1))
            molecule = match.group(2)
        else:
            coef = 1
            molecule = elet

        molecules[molecule] = ''
        result.append((coef, molecule))

    return result, molecules, rever


# print(molecule_sep('monohexosylceramide + 2 ferrocytochrome b5 + O2 + 2 H+'))
# print(molecule_sep('monohexosylceramide + ? + O2 + 2 H+'))
# print(molecule_sep('NAD+ + H+ |#116# 9% activity compared to cyclohexanone <197>| {r}'))
# print(molecule_sep('ATP + H2O {r}'))
# print(molecule_sep('ATP + H2O {ir}'))
# print(molecule_sep('2-(3-mol)test + mol2 + (S)-2-(3-mol)test'))

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
    RNX_data = []
    CMP_data = {}
    for element in data:
        reaction_SP = element['SP_data']
        i_symbol_egale = reaction_SP.find("=")
        substrates = molecule_sep(reaction_SP[:i_symbol_egale-1])
        produits = molecule_sep(reaction_SP[i_symbol_egale+2:])
        if produits == None or substrates == None:
            logging.warning('Exception')
        else:
            CMP_data.update(substrates[1])
            CMP_data.update(produits[1])
            # re-extrait les ID a partir de sub et prd
            #ou
            # mol sep retourne result et la list des ID pour faire le fichier CMP
            #enregistre avec une chaine vide et c'est plus tard que je mets le smile
            elets = element
            elets['reversibility'] = produits[2]
            elets['substrates'] = substrates[0]
            elets['products'] = produits[0]
            RNX_data.append(elets)

    with open(path + file_out, "w", encoding = 'utf8') as file:
        json.dump(RNX_data, file, indent = 2, ensure_ascii=False)
    file_out2 = 'out2.json'
    with open(path + file_out2, "w", encoding = 'utf8') as file:
        json.dump(CMP_data, file, indent = 2, ensure_ascii=False)

# modif_file(path, input_file, file_out)

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

def pubchem_cid_from_name(protein):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{protein}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID'][0]
    return None

# p = "ATP"
# pubchem_cid = pubchem_cid_from_name(p)

# url = 'https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl'
# client = Client(url)

# 
def mol_from_pubchem(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

# molfile_databis = mol_from_pubchem(pubchem_cid)
# print(molfile_databis)
# print(pubchem_cid)

def url_molfile_soap(email : str, password : str, prot_name):
    url = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    client = Client(url)
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    ligand_id = client.service.getLigandStructureIdByCompoundName(email, password, id=prot_name)
    if ligand_id:
        molfile_url = f'https://www.brenda-enzymes.org/molfile.php?LigandID={ligand_id}'
        return molfile_url
    else:
        return None

# url_molfile = url_molfile_soap("nolwenn.paris@inrae.fr",'brendamolfile', 'ATP')

def molfile_soap(url_molfile):
    response = requests.get(url_molfile)
    if response.status_code == 200:
        return response.text
    return None

# molfile_data = molfile_soap(url_molfile)

def molfile_to_smiles(molfile):
    try:
        mol = Chem.MolFromMolBlock(molfile)
        if mol is None:
            raise ValueError("nope premier etape de conversion en smile")
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        # print('le prbl est :', e)
        smiles = None

# smiles = molfile_to_smiles(molfile_data)
# print(smiles)

def file_mol_smile(path : str, input_file : str, file_out : str):
    with open(path + input_file, "r") as file:
        data = json.load(file)
    for protein in data.keys():
        try:
            url_molfile = url_molfile_soap("nolwenn.paris@inrae.fr",'brendamolfile', protein)
            molecule = molfile_soap(url_molfile)
            data[protein] = molfile_to_smiles(molecule)
        except Exception as e:
            print(e)
    with open(path + file_out, "w", encoding = 'utf8') as file:
        json.dump(data, file, indent = 2, ensure_ascii=False)

# file_mol_smile(path, 'out2.json', 'out3.json')
# =============================================================================
# PARTIE FUSION RXN / CMP
# =============================================================================
class RXN_CMP:
    def __init__(self, path, input_file, file_out1, file_out2):
        self.path = path
        self.input_file = input_file
        if file_out1:
            self.file_out_RXM = file_out1
        else:
            self.file_out = new_filename(input_file, 'RXN')
        if file_out2:
            self.file_out_CMP = file_out2
        else:
            self.file_out = new_filename(input_file, 'CMP')

    def get_path(self):
        return self.path
    def get_input_file(self):
        return self.input_file
    def get_file_out_RXN(self):
        return self.file_out_RXM
    def get_file_out_CMP(self):
        return self.file_out_CMP

    def run(self):
        modif_file(self.get_path(), self.get_input_file(), self.file_out_RXM())
        file_mol_smile(self.get_path(), self.file_out_RXM(), self.file_out_CMP())

# RXN_CMP(path, input_file, file_out).run()