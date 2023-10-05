#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 13:53:41 2023

@author: nparis
"""
from collections import defaultdict
from typing import Dict

path_data_brenda = '/home/nparis/brenda_enzyme/'
file_name_json = 'brenda_2023_1.json'
file_name_txt = 'brenda_2023_1.txt'

from brendapy import BrendaParser, BrendaProtein
from brendapy.console import console
from brendapy.log import get_logger
from brendapy.taxonomy import Taxonomy

BRENDA_PARSER = BrendaParser(str(path_data_brenda+file_name_txt))

# =============================================================================
# Element que nous voulons recuperer:
# EC number
# Km values
# Substrat
# Organisms
# Uniprot
# Commentary
# =============================================================================


# Recuperer la liste de tous les EC
list_EC_number = BRENDA_PARSER.keys()

# def parse_proteins_for_ec(ec: str = "1.1.1.1") -> Dict[int, BrendaProtein]:
#     """Parse the protein entries for a given EC number in BRENDA."""
#     proteins = BRENDA_PARSER.get_proteins(ec)
#     return proteins


# for ec_nb in list_EC_number:
#     parse_proteins_for_ec(ec_nb)

results=list()

p = BRENDA_PARSER.get_proteins("1.1.1.1")
for k, v in p.items():
    # if v.uniprot:
    results.append({"protein_id": k,"uniprot": p.uniprot,
                    "organism": p.organism,"KM": p.KM})

console.rule()
console.print(results)
console.rule()









































