#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Explication
"""
from argparse import ArgumentParser
from typing import Callable
# from os import getcwd as os_getcwd

# from brs_utils import add_logger_args

DEFAULT_ARGS = {
    'list_parameters' : ["ec", "uniprot", "organism", "substrate", 'comment', 'KM', 'TN', 'value'],
    'path_file_databrenda': '/home/nparis/brenda_enzyme/',
}


def build_args_parser(
        program : str,
        description : str,
        m_add_args: Callable = None) -> ArgumentParser:

    parser = ArgumentParser(
        program,
        description,
    )

    parser = add_arguments(parser)

    return parser


def add_arguments(parser : ArgumentParser) -> ArgumentParser:
    '''
    Build and return the ArgumentParser

    parser.add_argument(
        '',
        type= ,
        help='',
    )
    '''
    #Optional arg
    parser.add_argument(
        '--list_parameters',
        # action='append',
        nargs='+',
        default=DEFAULT_ARGS['list_parameters'],
        help='List of elements',
        # required=True
    )

    parser.add_argument(
        '--path_file_databrenda',
        type= str,
        default=DEFAULT_ARGS['path_file_databrenda'],
        help='Path to the Brenda file',
    )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debug mode'
        )
    return parser