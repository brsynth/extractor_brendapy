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


# path = '/home/nparis/brenda_enzyme/'
# input_file = 'test1.json'
# file_out = 'out1.json'

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
    tuple
        A tuple containing:
        - list of tuple: A list of tuples, where each tuple contains an integer coefficient 
          and a string molecule.
        - dict: A dictionary with molecule names as keys.
        - str: A string indicating whether the reaction is reversible ('r'), 
          irreversible ('ir'), or unknown ('?').
    """
    # Preprocessing: remove content within pipes and strip whitespace
    elements = re.sub(r'\|.*?\|', '', elements).strip()
    # Split the elements string by ' + ' to separate individual components
    eq = elements.split(' + ')
    molecules = {}
    result = []
    rever = '?'

    for elet in eq:
        elet = elet.strip()
        
        if elet == '?':
            return None

        # Determine the type of reaction
        if '{r}' in elet:
            rever = 'r'
            elet = elet.replace('{r}', '').strip()
        elif '{ir}' in elet:
            rever = 'ir'
            elet = elet.replace('{ir}', '').strip()

        # Match the coefficient and molecule name
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


def modif_file(path : str, input_file : str, file_out : str):
    """
    Modifies a JSON file by parsing chemical reactions and adding substrate 
    and product information.
    
    Parameters
    ----------
    path : str
        The directory path where the input and output files are located.
    input_file : str
        The name of the input JSON file containing reaction data.
    file_out : str
        The name of the output JSON file where modified data will be saved.
    
    Returns
    -------
    None
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
        # Parse substrates and products
        substrates = molecule_sep(reaction_SP[:i_symbol_egale-1])
        produits = molecule_sep(reaction_SP[i_symbol_egale+2:])
        if produits == None or substrates == None:
            logging.warning('Exception')
        else:
            # Update CMP_data with substrates and products
            CMP_data.update(substrates[1])
            CMP_data.update(produits[1])

            # Update reaction element with parsed data
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
    """
    Retrieve the PubChem Compound Identifier (CID) for a given protein name.

    This function queries the PubChem API to obtain the CID for a compound by its name.
    If the query is successful and a CID is found, the first CID from the results is returned.

    Parameters:
    -----------
    protein : str
        The name of the protein or compound to search for.

    Returns:
    --------
    int or None
        The CID of the compound if found, otherwise None.

    Raises:
    -------
    requests.exceptions.RequestException
        If there is an issue with the GET request.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{protein}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID'][0]
    return None


def mol_from_pubchem(cid):
    """
    Retrieve the molecular structure file (in SDF format) from PubChem for a given CID.

    This function queries the PubChem API to obtain the molecular structure file (SDF) for
    a compound by its CID. If the query is successful, the SDF content is returned.

    Parameters:
    -----------
    cid : int
        The PubChem Compound Identifier (CID) of the compound.

    Returns:
    --------
    str or None
        The molecular structure file content in SDF format if the request is successful,
        otherwise None.

    Raises:
    -------
    requests.exceptions.RequestException
        If there is an issue with the GET request.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

# molfile_databis = mol_from_pubchem(pubchem_cid)
# print(molfile_databis)
# print(pubchem_cid)

def url_molfile_soap(email : str, password : str, prot_name):
    """
    Retrieve the URL of the molecular file from the BRENDA Enzyme Database using SOAP.

    This function uses the BRENDA SOAP API to retrieve the ligand structure ID for a
    given protein name and constructs the URL for the molecular file.

    Parameters:
    -----------
    email : str
        The email address used for authentication with the BRENDA SOAP API.
    password : str
        The password used for authentication, which will be hashed before use.
    prot_name : str
        The name of the protein or ligand to search for.

    Returns:
    --------
    str or None
        The URL of the molecular file if the ligand ID is found, otherwise None.

    Raises:
    -------
    zeep.exceptions.Error
        If there is an issue with the SOAP request.
    """
    url = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    client = Client(url)
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    ligand_id = client.service.getLigandStructureIdByCompoundName(email, password, id=prot_name)
    if ligand_id:
        molfile_url = f'https://www.brenda-enzymes.org/molfile.php?LigandID={ligand_id}'
        return molfile_url
    else:
        return None


def molfile_soap(url_molfile):
    """
    Retrieve molecular file data from a given URL.

    This function sends a GET request to the specified URL to fetch the molecular file data.
    If the request is successful (HTTP status code 200), the content of the response is returned.
    Otherwise, None is returned.

    Parameters:
    -----------
    url_molfile : str
        The URL from which to retrieve the molecular file data.

    Returns:
    --------
    str or None
        The content of the molecular file if the request is successful, otherwise None.

    Raises:
    -------
    requests.exceptions.RequestException
        If there is an issue with the GET request.
    """
    response = requests.get(url_molfile)
    if response.status_code == 200:
        return response.text
    return None


def molfile_to_smiles(molfile):
    """
    Convert a molecular file format to SMILES notation.

    This function takes a molecular file in MOL format and converts it to SMILES notation using RDKit.
    If the conversion is successful, the SMILES string is returned. In case of an error during
    conversion, a warning is logged and None is returned.

    Parameters:
    -----------
    molfile : str
        The molecular file content in MOL format.

    Returns:
    --------
    str or None
        The SMILES string if the conversion is successful, otherwise None.

    Raises:
    -------
    ValueError
        If there is an error during the conversion process.
    """
    try:
        mol = Chem.MolFromMolBlock(molfile)
        if mol is None:
            raise ValueError("Error SMILE")
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        logging.warning(f'Exception : {e}')
        smiles = None

# smiles = molfile_to_smiles(molfile_data)

def file_mol_smile(path: str, input_file: str, file_out: str, mail : str, mdp :str, size = 100):
    """
    Process a JSON file containing protein data to convert molecular structures to SMILES notation.

    This function reads a JSON file, processes each protein entry by converting molecular
    structure data to SMILES notation, and writes the updated data to an output file.
    It also performs periodic temporary backups of the data during processing.

    Parameters:
    -----------
    path : str
        The directory path where the input and output files are located.
    input_file : str
        The name of the input JSON file containing protein data.
    file_out : str
        The name of the output JSON file where the processed data will be saved.
    mail : str
        Email address used for authentication with the external service.
    mdp : str
        Password used for authentication with the external service.
    size : int, optional
        The number of proteins to process before creating a temporary backup (default is 100).

    Returns:
    --------
    None

    Raises:
    -------
    Exception
        If an error occurs during the processing of a protein, the error is printed and
        the function continues with the next protein.

    Notes:
    ------
    - This function depends on the existence of `url_molfile_soap`, `molfile_soap`,
      and `molfile_to_smiles` functions, which are expected to handle external service
      communication and conversion logic.
    - Temporary backup files are created with the prefix 'temp_' followed by the output
      file name and saved in the specified path.
    """
    with open(path + input_file, "r") as file:
        data = json.load(file)
    
    count = 0
    
    for protein in data.keys():
        try:
            url_molfile = url_molfile_soap(mail, mdp, protein)
            molecule = molfile_soap(url_molfile)
            data[protein] = molfile_to_smiles(molecule)
        except Exception as e:
            print(e)
        
        count += 1
        
        if count % size == 0:
            temp_file_out = f"{path}temp_{file_out}"
            with open(temp_file_out, "w", encoding='utf8') as temp_file:
                json.dump(data, temp_file, indent=2, ensure_ascii=False)
            print(f"temporary backup {count}")

    with open(path + file_out, "w", encoding='utf8') as file:
        json.dump(data, file, indent=2, ensure_ascii=False)
    print('FINITO')

# file_mol_smile(path, 'out2.json', 'out3.json', "nolwenn.paris@inrae.fr", 'brendamolfile', size = 200)


# =============================================================================
# PARTIE FUSION RXN / CMP
# =============================================================================
class RXN_CMP:
    def __init__(self, path, input_file, file_out1, file_out2, mail, mdp):
        self.path = path
        self.input_file = input_file
        self.mail = mail
        self.mdp = mdp
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
    def get_mail(self):
        return self.mail
    def get_mdp(self):
        return self.mdp

    def run(self):
        modif_file(self.get_path(), self.get_input_file(), self.file_out_RXM())
        file_mol_smile(self.get_path(), self.get_file_out_RXM(),
                       self.get_file_out_CMP(), self.get_mail(), self.get_mdp())

