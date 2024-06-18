#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:52:30 2024

@author: nparis
"""
#conda install -c conda-forge rdkit
#pubchempy
#SOAP, RESTsolde

from rdkit import Chem

# from io import StringIO

from zeep import Client


# urlpath = 'https://www.brenda-enzymes.org/soap'
url = 'https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl'
client = Client(url)


# r_SOAP = client.service.get_molfile(enzyme_id='ATP') #1.1.1.1
# molfile_data = r_SOAP.molfile
# smiles = molfile_to_smiles(molfile_data)
# print(smiles)

import requests

def get_pubchem_cid_from_name(protein):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{protein}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID'][0]
    return None

p = "ATP"
pubchem_cid = get_pubchem_cid_from_name(p)

def get_molfile_from_pubchem(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

if pubchem_cid:
    molfile_data = get_molfile_from_pubchem(pubchem_cid)
    if molfile_data:
        print(molfile_data)
    else:
        print("Erreur avec pubmed")
else:
    print('nope')


def molfile_to_smiles2(molfile):
    try:
        # print(molfile)
        mol = Chem.MolFromMolBlock(molfile)
        print('mol', mol)
        if mol is None:
            raise ValueError("nope premier etape de conversion en smile")
        
        smiles = Chem.MolToSmiles(mol)
        
        return smiles
    
    except Exception as e:
        print('e', e)

smiles = molfile_to_smiles2(molfile_data)
print(smiles)