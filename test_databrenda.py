#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:50:57 2023

@author: nparis
"""
import unittest
import testbrendapy

# modifier les test pour le faire sur les plus petit fonction
class TestDataBrenda(unittest.TestCase):

    def test_good_nb_data(self):
        setdataprot = testbrendapy.data_brenda(['1.1.1.1'], 'KM')
        self.assertEqual(len(setdataprot), 299)

    def test_nb_organisms(self):
        setdataprot = testbrendapy.data_brenda(['1.1.1.1'], 'KM')
        count = 0
        for d_prot in setdataprot:
            if d_prot['organism'] == 'Homo sapiens':
                count += 1
        self.assertEqual(count,48)

    def test_name_new_file_created(self):
        self.assertEqual(testbrendapy.name_new_file_created('KKM'), 'setbrenda_KKM.json')

if __name__ == '__main__':
    unittest.main()
