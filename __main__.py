#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""

import sys
from .args import build_args_parser
from .datasetbrenda import DataSetBrenda


#DataSetBrenda(list_p, '/home/nparis/brenda_enzyme/').run()
#Mettre en premier la liste puis le chemin


def main():
    """
    Main function
    """
    try:
        opts, arg = sys.getopt.getopt(sys.argv, "hg:d", ["help", "grammar="])
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

    parser = build_args_parser(
        program='',
        description=''
    )

    args = parser.parse_args(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main(sys.argv[1:])