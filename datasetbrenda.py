#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:10:09 2023

@author: nparis
"""
import testbrendapy
# from brendapy import BrendaParser, BrendaProtein

class DataSetBrenda:
    def __init__(self, list_paramater : list, path_data_brenda : str,
                 list_ec : list = []):
        self.d_parameter_setting = testbrendapy.parameter_sorting(list_paramater)
        #noms du fichier avec la date et heure de creation ...
        self.path_set_brend = path_data_brenda + testbrendapy.name_new_file_created()

        if list_ec:
            self.list_ec = list_ec
        else:
            self.list_ec = testbrendapy.list_all_ec_in_data()

    def get_cinetique_parameter(self):
        return self.d_parameter_setting
    def get_path_set_brend(self):
        return self.path_set_brend
    def get_list_ec(self):
        return self.list_ec

    #@property -> enlever () apres le run dans l'appel de fonction
    def run(self):
        testbrendapy.create_file_json(self.get_path_set_brend(),
                         testbrendapy.data_brenda(self.get_list_ec(),
                                     self.get_cinetique_parameter()))