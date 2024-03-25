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

    parser = build_args_parser(
        program = 'extractor_brendapy',
        description = '',
        m_add_args = add_arguments
    )
    args = parser.parse_args()


    # laisser le run tant qu'il n'y a pas @property
    DataSetBrenda(args.list_parameters, args.path_file_databrenda, args.list_ec).run()
    # print(args)

if __name__ == "__main__":
    main()