#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""
import unittest
import os
import json
from tempfile import TemporaryDirectory
import RMX_file_creator as RMX

class TestRMXFile(unittest.TestCase):
    def test1_molecule_sep(self):
        r = RMX.molecule_sep('monohexosylceramide + 2 ferrocytochrome b5 + O2 + 2 H+')
        desired_r = [(1, 'monohexosylceramide'), (2, 'ferrocytochrome b5'),
                     (1, 'O2'), (2, 'H')]
        self.assertListEqual(r, desired_r)

    def test2_molecule_sep(self):
        r = RMX.molecule_sep('ATP + H2O')
        desired_r = [(1, 'ATP'), (1, 'H2O')]
        self.assertListEqual(r, desired_r)
    
    def test_add_rxn_filename(self):
        self.assertEqual('RMX_test1.json', RMX.add_rmx_filename('test1.json'))

if __name__ == '__main__':
    unittest.main()