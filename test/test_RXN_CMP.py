#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""
import unittest
import os
import json
from tempfile import TemporaryDirectory
import RXN_CMP_file_creator as RNX_CMP

class TestRXNFile(unittest.TestCase):
    def test1_molecule_sep(self):
        r = RNX_CMP.molecule_sep('monohexosylceramide + 2 ferrocytochrome b5 + O2 + 2 H+')
        desired_r = ([(1, 'monohexosylceramide'), (2, 'ferrocytochrome b5'),
                     (1, 'O2'), (2, 'H+')] , {'monohexosylceramide': '',
                                             'ferrocytochrome b5': '', 'O2': '',
                                             'H+': ''}, '?')
        self.assertEqual(r, desired_r)

    def test2_molecule_sep(self):
        r = RNX_CMP.molecule_sep('ATP + H2O {r}')
        desired_r = ([(1, 'ATP'), (1, 'H2O')], {'ATP': '', 'H2O': ''}, 'r')
        self.assertEqual(r, desired_r)
    
    def test3_molecule_sep(self):
        r = RNX_CMP.molecule_sep('2-(3-mol)test + mol2 + (S)-2-(3-m)test')
        desired_r = ([(1, '2-(3-mol)test'), (1, 'mol2'), (1, '(S)-2-(3-m)test')],
                     {'2-(3-mol)test': '', 'mol2': '', '(S)-2-(3-m)test': ''}, '?')
        self.assertEqual(r, desired_r)

    def test4_molecule_sep(self):
        r = RNX_CMP.molecule_sep('monohexosylceramide + ? + O2 + 2 H+')
        desired_r = (None)
        self.assertEqual(r, desired_r)
    
    def test5_molecule_sep(self):
        r = RNX_CMP.molecule_sep('NAD+ + H+ |#116# 9% activity compared to cyclohexanone <197>| {ir}')
        desired_r = ([(1, 'NAD+'), (1, 'H+')], {'NAD+': '', 'H+': ''}, 'ir')
        self.assertEqual(r, desired_r)
        
    # def test_modif_file(self):
    #     # Create temporary files for testing
    #     test_path = './test_data/'
    #     input_file = 'input.json'
    #     output_file = 'output.json'
    #     cmp_file = 'cmp_data.json'

    #     if not os.path.exists(test_path):
    #         os.makedirs(test_path)

    #     input_data = [{'SP_data': '2 H2O + 1 O2 = 2 H2O2'}]

    #     with open(test_path + input_file, 'w') as f:
    #         json.dump(input_data, f, indent=2)

    #     # Run the function
    #     RNX_CMP.modif_file(test_path, input_file, output_file)

    #     # Check output files
    #     with open(test_path + output_file, 'r') as f:
    #         output_data = json.load(f)

    #     with open(test_path + cmp_file, 'r') as f:
    #         cmp_data = json.load(f)

    #     self.assertEqual(len(output_data), 2)
    #     self.assertIn('reversibility', output_data[0])
    #     self.assertIn('substrates', output_data[0])
    #     self.assertIn('products', output_data[0])
    #     self.assertGreater(len(cmp_data), 0)

    #     # Cleanup
    #     os.remove(test_path + input_file)
    #     os.remove(test_path + output_file)
    #     os.remove(test_path + cmp_file)
    #     os.rmdir(test_path)
    
    def test1_new_filename(self):
        self.assertEqual('RXN_test1.json', RNX_CMP.new_filename('test1.json', 'RXN'))

    def test2_new_filename(self):
        self.assertEqual('CPM_test1.json', RNX_CMP.new_filename('test1.json', 'CPM'))

# class TestCMPFile(unittest.TestCase):
#     pass

if __name__ == '__main__':
    unittest.main()

