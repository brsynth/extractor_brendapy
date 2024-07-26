# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from typing import Callable


modif_file(self.get_path(), self.get_input_file(), self.file_out_RXM())
file_mol_smile(self.get_path(), self.file_out_RXM(),
               self.file_out_CMP(), self.mail(), self.mdp())

DEFAULT_ARGS = {
    'path' : '/home/nparis/brenda_enzyme/',
    'input_file' : 'file1',
    'file_out_RXM' : 'fileRXM',
    'file_out_CMP' : 'fileCMP',
    'mail' : '',
    'mdp' : ''
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
    parser.add_argument(
        '--path_file',
        type= str,
        default=DEFAULT_ARGS['path'],
        help='Path to the Brenda file',
    )

    parser.add_argument(
        '--input_file',
        type= str,
        default=DEFAULT_ARGS['input_file'],
        help='file all parameter',
    )
    parser.add_argument(
        '--file_RXM',
        type= str,
        default=DEFAULT_ARGS['file_out_RXM'],
        help='file RXM parameter',
    )
    parser.add_argument(
        '--file_CMP',
        type= str,
        default=DEFAULT_ARGS['file_out_CMP'],
        help='file smile',
    )
    
    parser.add_argument(
        '--mail',
        type= str,
        default=DEFAULT_ARGS['mail'],
        help='mail for brenda',
    )
    parser.add_argument(
        '--mdp',
        type= str,
        default=DEFAULT_ARGS['mdp'],
        help='password for brenda',
    )

    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debug mode'
        )
    return parser
