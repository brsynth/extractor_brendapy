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


if __name__ == '__main__':
    unittest.main()
