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
                     (1, 'O2'), (2, 'H')] , {'monohexosylceramide': '',
                                             'ferrocytochrome b5': '', 'O2': '',
                                             'H': ''})
        self.assertEqual(r, desired_r)

    def test2_molecule_sep(self):
        r = RNX_CMP.molecule_sep('ATP + H2O')
        desired_r = ([(1, 'ATP'), (1, 'H2O')], {'ATP': '', 'H2O': ''})
        self.assertEqual(r, desired_r)
    
    def test3_molecule_sep(self):
        r = RNX_CMP.molecule_sep('2-(3-mol)test + mol2 + (S)-2-(3-mol)test')
        desired_r = [(1, '2-(3-mol)test'), (1, 'mol2'), (1, '(S)-2-(3-mol)test')]
        self.assertEqual(r, desired_r)

    def test1_new_filename(self):
        self.assertEqual('RXN_test1.json', RNX_CMP.new_filename('test1.json', 'RXN'))

    def test2_new_filename(self):
        self.assertEqual('CPM_test1.json', RNX_CMP.new_filename('test1.json', 'CPM'))

# class TestRXNFile(unittest.TestCase):
#     pass

if __name__ == '__main__':
    unittest.main()