#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:50:57 2023

@author: nparis
"""
import unittest
import testbrendapy

class TestDataBrenda(unittest.TestCase):

    def test_file_path_request(self):
        result = testbrendapy.file_path_request('/home/nparis/brenda_enzyme/',
                                                'brenda_2023_1.txt')
        self.assertEqual(result, '/home/nparis/brenda_enzyme/brenda_2023_1.txt')

    def test2_file_path_request(self):
        pass

    def test_list_all_ec_in_data(self):
        pass
    def test_is_parameter_values(self):
        pass
    def test_find_shared_substrate(self):
        pass
    def test_d_comment_each_kinetic(self):
        pass

    def test_find_keys_with_similar_values(self):
        d_test = {'TN': {'16': 'pH 7.0, 25°C, mutant N107D <17>'},
                  'KM': {'20': 'pH 7.0, 25°C, mutant N107D <17>',
                         '21': 'pH 7.0, 25°C, mutant N107L <17>'}}
        result = [{'TN': '16', 'KM': '20'}]
        self.assertEqual(testbrendapy.find_keys_with_similar_values(d_test), result)

    def test2_find_keys_with_similar_values(self):
        d_test = {'TN': {'16': 'mutant N107D'},
                  'KM': {'20': 'mutant N107D',
                         '21': 'mutant N107L'},
                  'KKM': {'2': 'mutant N107D',
                         '1': 'mutant N107L'}}
        result = [{'TN': '16', 'KM': '20', 'KKM': '2'}, {'KM': '21', 'KKM': '1'}]
        self.assertEqual(testbrendapy.find_keys_with_similar_values(d_test), result)

    def test_create_subdict_json(self):
        pass

    def test_commun_lists(self):
        l1 = ['un', 'deux', 'trois']
        l2 = ['un', 'quatre']
        self.assertEqual(testbrendapy.commun_lists(l1,l2), ['un'])

    def test2_commun_lists(self):
        l1 = ['deux', 'trois']
        l2 = ['un', 'quatre']
        self.assertEqual(testbrendapy.commun_lists(l1,l2), [])

    def test_parameter_sorting(self):
        l1 = ['ec', 'uniprot', 'value', 'units', 'KM']
        result = {'p_str': ['uniprot', 'ec'],
                  'p_list_dict': ['KM'],
                  'key_p_list_dict': ['value', 'units'],
                  'p_set': []}
        self.assertDictEqual(testbrendapy.parameter_sorting(l1), result)

    def test2_parameter_sorting(self):
        l2 = ['ec', 'uniprot', 'value', 'units', 'KM', 'IC50', 'tissues']
        result = {'p_str': ['uniprot', 'ec'],
                  'p_list_dict': ['KM', 'IC50'],
                  'key_p_list_dict': ['value', 'units'],
                  'p_set': ['tissues']}
        self.assertDictEqual(testbrendapy.parameter_sorting(l2), result)

    # def test3_parameter_sorting(self):
    #     l3 = ['ec', 'uniprot', 'value', 'units', 'KM', 'tissues', 'intrus']
    #     result = 'Erreur'
    #     self.assertEqual(testbrendapy.parameter_sorting(l3), result)


if __name__ == '__main__':
    unittest.main()
