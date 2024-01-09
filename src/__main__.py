#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: noparis
"""

from .args import build_args_parser

def main():
    """
    Main function
    """
    parser = build_args_parser(
        program='',
        description=''
    )

    args = parser.parse_args()
