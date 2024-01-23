#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""

import sys
from .src.args import build_args_parser
from .src import DataSetBrenda

#Mettre en premier la liste puis le chemin
#Rajouter le lien ou les data genere doivnt etre depose 

def main():
    """
    Main function
    """
    try:
        opts, arg = sys.getopt.getopt(sys.argv, "hg:d", ["help"])
    except sys.getopt.GetoptError:
        #Si le script ne fonctionne pas
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            # usage()
            sys.exit()
        elif opt == '-d':
            global _debug
            _debug = 1

    # parser = build_args_parser()
    # args = parser.parse_args()

    parser = build_args_parser(
        program='',
        description=''
    )

    args = parser.parse_args(sys.argv[1], sys.argv[2])
    #args.list_parameter, args.path_data_brenda
    # DataSetBrenda(sys.argv[1], sys.argv[2]).run()

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        DataSetBrenda(sys.argv[1:]).run()
        #DataSetBrenda(sys.argv[1], sys.argv[2]).run()
    except TypeError:
        if not sys.argv[1:]:
            list_parameter : list = input('List of elements : ')
            path_data_brenda : str = input('Path to the Brenda file : ')
            DataSetBrenda(list_parameter, path_data_brenda).run()
    finally:
        print('The End')