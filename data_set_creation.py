#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:51:05 2023

@author: nparis
"""
from zeep import Client
import hashlib
import getpass
import collections

# =============================================================================
# Etape 1 : Connection à distance
# Etape 2 : Recuperation des donnees souhaitez
# Etape 3 : Vérifier que tous est bien obtenu (elimination de ceux ou il manque une donnée)
# Etape 4 : concationné et mettre au propres
# =============================================================================

def ask_password() -> str:
    """
    Demande le mdp pour ce connecté à distance à la BDD de brenda
    Attention le mdp n'est pas caché par des etoiles !!

    Returns : string qui contient le mdp
    -------
    str

    """
    return getpass.getpass('Password, please : ')

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256(ask_password().encode("utf-8")).hexdigest()
client = Client(wsdl)

# class Identifiant_connection(self):
#     pass
# Faire une demande pour identitfiant

parameters = ("nolwenn.paris@inrae.fr",password,"ecNumber*1.1.1.11",
              "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*",
              "id*", "organism*")
resultString = client.service.getSequence(*parameters)