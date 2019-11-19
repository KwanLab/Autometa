#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
function to retrieve orfs from provided assembly using prodigal
"""


import os
import subprocess


def get_orfs(assembly, nucls_out, prots_out, force=False, verbose=False):
    """Calls ORFs from provided input assembly
    Inputs:
        assembly - </path/to/assembly>
        nucls_out - </path/to/nucls.out>
        prots_out - </path/to/prots.out>
        force - True|False
    ReturnType: 3-tuple (boolean, str, str)
    Returns:
        On Success - (True, nucls_out, prots_out)
        On Failure - (False, None, None)
    """
    # OPTIMIZE: This should use parallel or multiprocessing
    # OPTIMIZE: Piped for parallel/multiprocessing rather than writing first to file
    if not os.path.exists(assembly):
        raise FileNotFoundError(f'{assembly} does not exists!')
    for fpath in [nucls_out, prots_out]:
        if os.path.exists(nucls_out) and not force:
            raise FileExistsError(f'{fpath} To overwrite use --force')
    cmd = [
        'prodigal',
        '-i',assembly,
        '-a',prots_out,
        '-d',nucls_out,
        '-p',
        'meta',
        '-m',
        '-q'
    ]
    if verbose:
        print(f'RunningProdigal: {" ".join(cmd)}')
    with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
        proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
    if proc.returncode:
        if verbose:
            print(f'ProdigalFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
        return (False, None, None)
    else:
        return (True, nucls_out, prots_out)

def main(args):
    success, nucls_out, prots_out = get_orfs(
        assembly=args.assembly,
        nucls_out=args.nucls_out,
        prots_out=args.prots_out,
        force=args.force,
        verbose=args.verbose)
    if success:
        if args.verbose:
            print(f'written:\nnucls fpath: {nucls_out}\nprots fpath: {prots_out}')
        return 0
    else:
        if args.verbose:
            print('Prodigal ORF calling failed')
        return 1

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Calls ORFs with provided input assembly')
    parser.add_argument('assembly', help='</path/to/assembly>')
    parser.add_argument('nucls_out', help='</path/to/nucls.out>')
    parser.add_argument('prots_out', help='</path/to/prots.out>')
    parser.add_argument('--force', help="force overwrite of ORFs out filepaths",
        default=True)
    parser.add_argument('--verbose', help="add verbosity", default=True)
    args = parser.parse_args()
    main(args)
