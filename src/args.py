#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Explication
"""
from argparse import (
    ArgumentParser,
    # BooleanOptionalAction
)
from os import getcwd as os_getcwd

from brs_utils import add_logger_args

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
    

    Parameters
    ----------
    parser : TYPE
        DESCRIPTION.

    Returns
    -------
    parser : TYPE
        DESCRIPTION.

    '''
    parser.add_argument(
        'list_parameters',
        type= list,
        help='list parameters for data set',
    )

    parser.add_argument(
        'path_file_databrenda',
        type= str,
        help='',
    )
    return parser