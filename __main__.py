#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""

import sys
from .src.args import add_arguments, build_args_parser
from .src import DataSetBrenda
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger


def main():
    """
    Main function
    """
    # try:
    #     opts, arg = sys.getopt.getopt(sys.argv, "hg:d", ["help"])
    # except sys.getopt.GetoptError:
    #     #Si le script ne fonctionne pas
    #     sys.exit(2)
    
    # for opt, arg in opts:
    #     if opt in ("-h", "--help"):
    #         # usage()
    #         sys.exit()
    #     elif opt == '-d':
    #         global _debug
    #         _debug = 1
    
    parser = build_args_parser(
        program = 'brenda_enz_code',
        description = '',
        m_add_args = add_arguments
    )
    args = parser.parse_args()
    # print(args.list_parameters)
    # vars(args)
    # laisser le run tant qu'il n'y a pas @property
    DataSetBrenda(args.list_parameters, args.path_file_databrenda).run()

if __name__ == "__main__":
    main()