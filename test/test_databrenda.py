#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:50:57 2023

@author: nparis
"""
import unittest
import testbrendapy
from collections import OrderedDict

class TestDataBrenda(unittest.TestCase):

    def test_file_path_request(self):
        result = testbrendapy.file_path_request('/home/nparis/brenda_enzyme/',
                                                'brenda_2023_1.txt')
        self.assertEqual(result, '/home/nparis/brenda_enzyme/brenda_2023_1.txt')

    def test_is_parameter_values(self):
        dict_test = OrderedDict([('test1', 5), ('test2', '1.1.1.103'),
                                 ('test3', 'Homo sapiens'), ('test4', None),
                                 ('test5', {1: {'t1': 'b1', 't2': 123}}),
                                 ('test6', set())])
        list_test = ['test1','test2', 'test3', 'test5']
        self.assertTrue(testbrendapy.is_parameter_values(list_test, dict_test))

    def test2_is_parameter_values(self):
        dict_test = OrderedDict([('test1', 5), ('test2', '1.1.1.103'),
                                 ('test3', 'Homo sapiens'), ('test4', None),
                                 ('test5', {1: {'t1': 'b1', 't2': 123}}),
                                 ('test6', set())])
        list_test = ['test2', 'test4']
        self.assertFalse(testbrendapy.is_parameter_values(list_test, dict_test))

    def test3_is_parameter_values(self):
        dict_test = OrderedDict([('test1', 5), ('test2', '1.1.1.103'),
                                 ('test3', 'Homo sapiens'), ('test4', None),
                                 ('test5', {1: {'t1': 'b1', 't2': 123}}),
                                 ('test6', set())])
        list_test = ['test1','test5', 'test6']
        self.assertFalse(testbrendapy.is_parameter_values(list_test, dict_test))

    def test_find_shared_substrate(self):
        parameter = 'TN'
        d_test = [{'comment': 'mutant enzyme M333E, at pH 7.5 and 37°C <46>',
                   'value': 46.7, 'substrate': 'L-threonine'},
                  {'comment': 'wild type enzyme, at pH 7.5 and 37°C <46>',
                   'value': 47.3, 'substrate': 'L-threonine'}]
        d_result = {'L-threonine': {'TN': [0, 1]}}
        self.assertDictEqual(testbrendapy.find_shared_substrate({}, d_test, parameter), d_result)

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
        result = {'p_str': ['ec','uniprot'],
                  'p_list_dict': ['KM'],
                  'key_p_list_dict': ['value', 'units'],
                  'p_set': []}
        self.assertDictEqual(testbrendapy.parameter_sorting(l1), result)

    def test2_parameter_sorting(self):
        l2 = ['ec', 'uniprot', 'value', 'units', 'KM', 'IC50', 'tissues']
        result = {'p_str': ['ec','uniprot'],
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
