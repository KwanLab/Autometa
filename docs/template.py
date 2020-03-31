#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.

Template Script for Autometa Modules

Template Description:
Concise sentence describing functionality of script

Template Format:
0. - Shebang python env definition
1. - Two lines following comment block script description
2. - One line separating different import types
3. - Two lines separating imports and logger
4. - Two lines separating logger and algorithm functions
5. - Main function
6. - if __name__ == '__main__' clause
7. - Argparse
8. - Logging aliased to logger in clause 6.
9. - Pass args to main
"""


import logging
# import <package>


# import <package> as <package-alias>

# from <library> import <moduleA>
# from <library> import <moduleB>

# from <autometa.common> import <module>
# from <autometa> import <module>
# from <autometa.module> import <function>


# Two lines to separate imports and beginning of script
logger = logging.getLogger(__name__)
# Two lines to separate function definitions and logger


def fname(arg):
    """Short summary.

    Parameters
    ----------
    arg : type
        Description of parameter `arg`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    pass

def main(args):
    logger.info(args.hello_world)
    # operations on args.positional
    # operations on args.optional

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser(description='Concise Functional Description of Script')
    parser.add_argument('positional',help='<help text of positional arg>')
    parser.add_argument('--optional',help='<help text of optional arg>')
    parser.add_argument(
        '--hello-world',
        help='<help text of hello world>',
        default='Hello World')
    args = parser.parse_args()
    main(args)
