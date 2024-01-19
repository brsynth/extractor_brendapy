#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Explication
"""
from argparse import (
    ArgumentParser,
    # BooleanOptionalAction
)
# from os import getcwd as os_getcwd

# from brs_utils import add_logger_args

DEFAULT_ARGS = {
    'list_parameters' : ["ec", "uniprot", "organism", "substrate", 'comment', 'KM', 'TN', 'value'],
    'path_file_databrenda': '',
}


def build_args_parser(
        program,
        description):

    parser = ArgumentParser(
        program,
        description,
    )

    parser = add_arguments(parser)

    return parser


def add_arguments(parser):
    '''
    Build and return the ArgumentParser

    parser.add_argument(
        '',
        type= ,
        help='',
    )
    '''
    parser.add_argument(
        'list_parameters',
        type= list,
        help='List of elements',
    )

    parser.add_argument(
        'path_file_databrenda',
        type= str,
        help='Path to the Brenda file',
    )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debug mode'
        )
    return parser