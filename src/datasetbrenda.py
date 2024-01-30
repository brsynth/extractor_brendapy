#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:10:09 2023

@author: nparis
"""
from .testbrendapy import parameter_sorting, name_new_file_created
from .testbrendapy import data_brenda, create_file_json
from brendapy import BrendaParser

DEFAULT_FILE_NAME = 'brenda_2023_1.txt'

class DataSetBrenda:
    def __init__(self, list_paramaters : list, path_file_databrenda : str,
                 list_ec : list = []):
        self.d_parameter_setting = parameter_sorting(list_paramaters)
        #noms du fichier avec la date et heure de creation ...
        self.path_set_brend = path_file_databrenda + name_new_file_created()
        self.path_filebrenda = path_file_databrenda + DEFAULT_FILE_NAME
        # self.Brendadata = BrendaParser(self.path_filebrenda)

        # if list_ec:
        self.list_ec = list_ec
        # else:
        #     self.list_ec = list_all_ec_in_data(self.Brendadata)


    def get_cinetique_parameter(self):
        return self.d_parameter_setting
    def get_path_set_brend(self):
        return self.path_set_brend

    def get_list_ec(self, Brendadata):
        if self.list_ec:
            return self.list_ec
        else:
            self.list_ec = Brendadata.keys()
            return self.list_ec
    
    #@property -> enlever () apres le run dans l'appel de fonction
    def run(self):
        Brendadata = BrendaParser(self.path_filebrenda)
        create_file_json(self.get_path_set_brend(),
                         data_brenda(Brendadata, self.get_list_ec(Brendadata),
                                     self.get_cinetique_parameter()))



# list_p = ["ec", "uniprot"]
# DataSetBrenda(list_p, '/home/nparis/brenda_enzyme/').run()