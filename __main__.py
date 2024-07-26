#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""

import sys
from .src.args import add_arguments, build_args_parser
# from .src import DataSetBrenda
# from .src2.RXN_CMP_file_creator import RXN_CMP
from .src2.args import add_arguments, build_args_parser
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger


def main():
    """
    Main function
    """

    parser = build_args_parser(
        program = 'extractor_brendapy',
        description = '',
        m_add_args = add_arguments
    )
    args = parser.parse_args()


    # laisser le run tant qu'il n'y a pas @property
    DataSetBrenda(args.list_parameters, args.path_file_databrenda,
                  args.list_ec, args.namefile).run()
    # print(args)
    
    # parser = build_args_parser(
    #     program = 'RXN_CMP',
    #     description = '',
    #     m_add_args = add_arguments
    # )
    # args = parser.parse_args()

    # RXN_CMP.modif_file(args.path, args.input_file, args.file_RXM).run()
    # RXN_CMP.file_mol_smile(args.path, args.file_RXM, args.file_CMP, args.mail, args.mdp).run()

if __name__ == "__main__":
    main()