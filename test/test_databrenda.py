#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 13:50:57 2023

@author: nparis
"""
import unittest
from src import testbrendapy
from collections import OrderedDict
from datetime import datetime
import os
import json
from tempfile import TemporaryDirectory
from unittest import mock

class TestDataBrenda(unittest.TestCase):

    def setUp(self):
        self.dict_proteins = OrderedDict([('protein_id', 1),
                     ('ec', '1.1.1.98'),
                     ('organism', 'Homo sapiens'),
                     ('taxonomy', 9606),
                     ('uniprot', None),
                     ('CL',
                      [{'substrate' : 'substrate',
                        'data': '3 different genes, HAOX1, HAOX2, HAOX3',
                        'refs': [2]}]),
                     ('ID', '1.1.1.98'),
                     ('LO', [{'data': 'peroxisome', 'refs': [2]}]),
                     ('RN', {'(R)-2-hydroxy-fatty-acid dehydrogenase'}),
                     ('RT', {'redox reaction', 'reduction', 'oxidation'}),
                     ('references',
                      {1: {'info': 'Levis, G.M.: 2-Hydroxy fatty acid oxidases of rat '
                                   'kidney. Biochem. Biophys. Res. Commun. (1970) 38, '
                                   '470-477.',
                           'pubmed': 5443694}}),
                     ('tissues', {'BTO:0000988', 'BTO:0000671', 'BTO:0000759'})])

    def test_name_new_file_created(self):
        date_time = datetime.now()
        formatagedate = date_time.strftime('-%Y-%m-%d-%H-%M-%S')
        result = 'setbrenda_' + formatagedate + '.json'
        self.assertEqual(testbrendapy.name_new_file_created(), result)

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

    def test4_is_parameter_values(self):
        dict_test = OrderedDict([('test1', 5), ('test2', '1.1.1.103'),
                                  ('test3', 'Homo sapiens'), ('test4', None),
                                  ('test5', {1: {'t1': 'b1', 't2': 123}}),
                                  ('test6', set())])
        list_test = ['test1','test5', 'test7']
        self.assertFalse(testbrendapy.is_parameter_values(list_test, dict_test))

    def test_find_shared_substrate(self):
        parameter = 'TN'
        d_test = [{'comment': 'mutant enzyme M333E, at pH 7.5 and 37°C <46>',
                    'value': 46.7, 'substrate': 'L-threonine'},
                  {'comment': 'wild type enzyme, at pH 7.5 and 37°C <46>',
                    'value': 47.3, 'substrate': 'L-threonine'}]
        d_result = {'L-threonine': {'TN': [0, 1]}}
        self.assertDictEqual(testbrendapy.find_shared_substrate({}, d_test, parameter), d_result)

    def test2_find_shared_substrate(self):
        d_temporaire = {'L-threonine': {'KM': [0]}}
        parameter = 'TN'
        d_test = [{'comment': 'mutant enzyme M333E, at pH 7.5 and 37°C <46>',
                   'value': 46.7, 'substrate': 'L-threonine'},
                  {'comment': 'wild type enzyme, at pH 7.5 and 37°C <46>',
                    'value': 47.3, 'substrate': 'L-threonine'}]
        d_result = {'L-threonine': {'KM': [0], 'TN': [0, 1]}}
        self.assertDictEqual(testbrendapy.find_shared_substrate(d_temporaire, d_test,
                                                                parameter), d_result)

    def test_d_comment_each_kinetic(self):
        d_index = {}
        d_i_substr = {'TN': [16], 'KM': [20, 21]}
        dict_proteins = {'TN': {16: {'comment': 'pH 7.0, 25°C, mutant N107D <17>'}},
                        'KM': {20: {'comment': 'pH 7.0, 25°C, mutant N107D <17>'},
                                21: {'comment': 'pH 7.0, 25°C, mutant N107L <17>'}}}
        result = {'TN': {'16': 'pH 7.0, 25°C, mutant N107D <17>'},
                  'KM': {'20': 'pH 7.0, 25°C, mutant N107D <17>',
                        '21': 'pH 7.0, 25°C, mutant N107L <17>'}}
        self.assertEqual(testbrendapy.d_comment_each_kinetic(d_index, d_i_substr,
                                                            dict_proteins), result)

    def test2_d_comment_each_kinetic(self):
        d_index = {}
        d_i_substr = {'KM': [20, 21]}
        dict_proteins = {'KM': {20: {'comment': 'comment1'}, 21: {}}}
        
        result = {'KM': {'20': 'comment1'}}
        self.assertEqual(testbrendapy.d_comment_each_kinetic(d_index, d_i_substr,
                                                              dict_proteins), result)

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
        d_temporaire = {}
        d_p_setting = {'p_str': ['param1', 'param2'],
                        'key_p_list_dict': ['key1', 'key2']}
        dict_proteins = {'param1': 'value1',
                          'param2': 'value2',
                          'kinetic_value': {'sub_d_brenda_value':
                                            {'key1': 'result1',
                                            'key2': 'result2'}}}
        p_kinetic = 'kinetic_value'
        i_sub_d_brenda = 'sub_d_brenda_value'

        test = testbrendapy.create_subdict_json(d_temporaire, d_p_setting,
                                                dict_proteins, i_sub_d_brenda,
                                                p_kinetic)
        result = {'param1': 'value1','param2': 'value2',
                  'kinetic_value_key1': 'result1','kinetic_value_key2': 'result2'}

        self.assertEqual(test, result)

    def test2_create_subdict_json(self):
        d_temporaire = {}
        d_p_setting = {'p_str': ['param1'], 'key_p_list_dict': ['key1']}
        dict_proteins = {'param1': 'value1'}
        i_sub_d_brenda = 'kinetic_value'
        p_kinetic = 'sub_d_brenda_value'
        
        result = {'param1': 'value1'}

        self.assertEqual(testbrendapy.create_subdict_json(d_temporaire, d_p_setting,
                                             dict_proteins, i_sub_d_brenda,
                                             p_kinetic), result)

    def test_create_file_json(self):
        with TemporaryDirectory() as temp_dir:
            # Repertoire temporaire
            json_path = os.path.join(temp_dir, "test.json")

            test_data = [{"1.1.1.1": "ec", "substrat": "name1"},
                          {"2.2.2.2": "ec", "substrat": "name2"},
                          {"3.3.3.3": "ec", "substrat": "name3"}]

            testbrendapy.create_file_json(json_path, test_data)
            self.assertTrue(os.path.exists(json_path))

            with open(json_path, "r", encoding='utf-8') as file:
                loaded_data = json.load(file)
            self.assertEqual(test_data, loaded_data)

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
        d_test = testbrendapy.parameter_sorting(l1)
        for key in result :
            self.assertListEqual(sorted(d_test[key]), sorted(result[key]))

    def test2_parameter_sorting(self):
        l2 = ['ec', 'uniprot', 'value', 'units', 'KM', 'IC50', 'tissues']
        result = {'p_str': ['ec','uniprot'],
                  'p_list_dict': ['KM', 'IC50'],
                  'key_p_list_dict': ['value', 'units'],
                  'p_set': ['tissues']}
        d_test = testbrendapy.parameter_sorting(l2)
        for key in result :
            self.assertListEqual(sorted(d_test[key]), sorted(result[key]))

    def test_check_parameter_values(self):
        d_p_setting = {'p_str': ['test1'], 'p_list_dict': ['test2', 'test3']}
        d1 = {'test1': '1', 'test2': 'v1', 'test3': 'v2'}
        result = testbrendapy.check_parameter_values(d_p_setting, d1)
        self.assertTrue(result)

    def test2_check_parameter_values(self):
        d_p_setting = {'p_str': ['test1'], 'p_list_dict': ['test2', 'test3']}
        d2 = {'test1': '1', 'test2': 'v1'}
        result = testbrendapy.check_parameter_values(d_p_setting, d2)
        self.assertFalse(result)

    def test_find_shared_substrate_index(self):
        para_list = ['test1']
        d1 = {'test1': {'k1': 1, 'k2': 2}}
        result = testbrendapy.find_shared_substrate_index(para_list, d1)
        self.assertEqual(result,  {})

    def test2_find_shared_substrate_index(self):
        para_list = ['test1', 'test2']
        d2 = {'test1': [{'substrate': 1}, {'k2': 2}],
              'test2': [{'k2': 3}, {'substrate': 4}]}
        result = testbrendapy.find_shared_substrate_index(para_list, d2)
        self.assertEqual(result, {'1': {'test1': [0]}, '4': {'test2': [1]}})

    def test3_find_shared_substrate_index(self):
        para_list = ['test1', 'test2', 'test3']
        d3 = {'test1': [{'substrate': 1}, {'k2': 2}],
              'test2': [{'k2': 3}, {'substrate': 4}],
              'test3' : [{'k3': 3}, {'k1': 5}, {'k2': 6}, {'substrate': 4}]}
        result = testbrendapy.find_shared_substrate_index(para_list, d3)
        self.assertEqual(result, {'1': {'test1': [0]}, '4': {'test2': [1], 'test3': [3]}})

    def test_pre_subdict_from_couple(self):
        d_p_setting = {'p_list_dict': []}
        d1 = {}
        couple = {}
        result = testbrendapy.pre_subdict_from_couple(d_p_setting, d1, couple)
        self.assertEqual(result, {})

    def test2_pre_subdict_from_couple(self):
        d_p_setting = {'p_str': ['param1'],
                        'p_list_dict': ['test1', 'test2'],
                        'key_p_list_dict' : ['key1']}
        dict_proteins = {'param1': 'value1',
                         'test1': [{'key1': '1', 'k2': '2'}, {'key2': '1', 'k2': '2'}],
                         'test2': [{'key3': '1', 'k2': '2'}, {'key1': '1', 'k2': '2'}]}
        couple = {'test1': 0}
        result = testbrendapy.pre_subdict_from_couple(d_p_setting, dict_proteins, couple)
        self.assertEqual(result, {'param1': 'value1', 'test1_key1': '1'})

    def test3_pre_subdict_from_couple(self):
        d_p_setting = {'p_str': ['param1'],
                        'p_list_dict': ['test1', 'test2'],
                        'key_p_list_dict' : ['key1']}
        dict_proteins = {'param1': 'value1',
                         'test1': [{'key1': '1', 'k2': '2'}, {'key2': '1', 'k2': '2'}],
                         'test2': [{'key3': '1', 'k2': '2'}, {'key1': '4', 'k2': '2'}]}
        couple = {'test1': 0, 'test2': 1}
        result = testbrendapy.pre_subdict_from_couple(d_p_setting, dict_proteins, couple)
        # 'test2_key1': '1'
        self.assertEqual(result, {'param1': 'value1', 'test1_key1': '1', 'test2_key1': '4'})

    def test1_data_brenda(self):
        list_ec = ['1.1.1.98']
        d_p_setting = {'p_str': ["ec", "organism"],
                        'p_list_dict': ['CL'],
                        'key_p_list_dict' : ['data']
                        }
        # Librairie Mock
        #

        # a = mock.MagicMock(name = 'dict_proteins')
        # a.data.return_value = self.dict_proteins
        # print(a.data())

        mock_brenda_parser2 = mock.MagicMock(name = 'dict_proteins')
        mock_brenda_parser = mock.Mock()
        # mock_brenda_parser2.data = Mock(return_value = self.dict_proteins)
        # mock_brenda_parser = MagicMock(get_proteins = {'data' : mock_brenda_parser2})
        mock_brenda_parser.get_proteins = mock.Mock(return_value = {'data' : mock_brenda_parser2})
        mock_brenda_parser2.data = self.dict_proteins
        results = testbrendapy.data_brenda(mock_brenda_parser, list_ec, d_p_setting)
        result_expected = [{'ec': '1.1.1.98',
                            'organism': 'Homo sapiens',
                            'CL_data': '3 different genes, HAOX1, HAOX2, HAOX3'}]
        # print(mock_brenda_parser.get_proteins())
        # print('_____________')
        # print('r', results)
        self.assertEqual(results, result_expected)
        # rbis = testbrendapy.data_brenda(Brendadata, ['1.1.1.98'], d_p_setting)
        # print('...', rbis)

    # def test2(self):
    #     from brendapy import BrendaParser
    #     Brendadata = BrendaParser('/home/nparis/brenda_enzyme/brenda_2023_1.txt')
    #     d_p_setting = {'p_str': ["ec", "organism", 'uniprot'],
    #                     'p_list_dict': ['KM'],
    #                     'key_p_list_dict' : ['value', 'substrate', 'comment']
    #                     }
    #     rbis = testbrendapy.data_brenda(Brendadata, ['1.1.1.1'], d_p_setting)
    #     print('...', rbis)


if __name__ == '__main__':
    unittest.main()