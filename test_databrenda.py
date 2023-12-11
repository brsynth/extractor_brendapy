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

    def test_find__keys_with_similar_values(self):
        d_test = {'TN': {'16': 'pH 7.0, 25°C, mutant N107D <17>'},
                  'KM': {'20': 'pH 7.0, 25°C, mutant N107D <17>',
                         '21': 'pH 7.0, 25°C, mutant N107L <17>'}}
        result = [{'TN': '16', 'KM': '20'}]
        self.assertEqual(testbrendapy.find__keys_with_similar_values(d_test), result)

    def test_create_subdict_json(self):
        pass
    def test_create_file_json(self):
        pass
    def test_commun_lists(self):
        pass
    def test_parameter_sorting(self):
        pass


if __name__ == '__main__':
    unittest.main()
