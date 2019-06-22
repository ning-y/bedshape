#! /usr/bin/env python3

import argparse, logging, re, sys

logger = logging.getLogger('bedshape')

def get_parser():
    parser = argparse.ArgumentParser('bedshape')
    subparser_adder = parser.add_subparsers()
    parser.set_defaults(which_subcommand=None)

    alias_parser = subparser_adder.add_parser(
            'alias', description='Manage aliases for filenames.')
    config_alias(alias_parser)
    profile_parser = subparser_adder.add_parser(
            'profile', description='Calculate reactivity profiles.')
    config_profile(profile_parser)

    return (parser, alias_parser, profile_parser)

def get_root_parser():
    return get_parser()[0]

def config_alias(parser):
    parser.set_defaults(which_subcommand='alias')
    subparser_adder = parser.add_subparsers()

    list_parser = subparser_adder.add_parser('list',
            description='List all available aliases.')
    list_parser.set_defaults(which_subcommand='list')

    add_parser = subparser_adder.add_parser('add',
            description='Add a new alias.')
    add_parser.set_defaults(which_subcommand='add')
    add_parser.add_argument('name', type=str,
            help='Name of alias to add')
    add_parser.add_argument('--reference', '-r', required=True, type=str,
            help='Absolute path to reference for this alias.')
    add_parser.add_argument('--modified', '-m', required=True, type=str,
            help='Absolute path to the modified alignment for this alias.')
    add_parser.add_argument('--unmodified', '--untreated', '-u', type=str,
            help='Absolute path to the unmodified alignment for this alias.')
    add_parser.add_argument('--denatured', '-d', type=str,
            help='Absolute path to the denatured alignment for this alias.')

    rm_parser = subparser_adder.add_parser('rm',
            description='Remove aliases.')
    rm_parser.set_defaults(which_subcommand='rm')
    rm_parser.add_argument('name', type=str,
            help='Name of alias to remove')

def config_profile(parser):
    parser.set_defaults(which_subcommand='profile')

    parser.add_argument('--alias', '-a', type=str,
            help='Alias to use. This, or a combination of --reference and '
                '--modified, is required.')
    parser.add_argument('--reference', '-r', type=str,
            help='Reference to use. A combination of this and --modified, or '
                '--alias, is required.')
    parser.add_argument('--modified', '-m', type=str,
            help='Modified alignment to use. A combination of this and '
                '--reference, or --alias, is required.')
    parser.add_argument('--unmodified', '-u', type=str,
            help='Unmodified alignment to use.')
    parser.add_argument('--denatured', '-d', type=str,
            help='Denatured alignment to use.')

    parser.add_argument('--region', '-rg', type=str,
            help='Region to use. <rname>:<start>-<stop>. Commas are allowed '
                'in <start> and <stop>. This, or --bed is required.')
    parser.add_argument('--bed', '-b', type=str,
            help='BED file containing regions to use. This, or --region '
                'is required.')

    parser.add_argument('--min-depth', type=int, default=0,
            help="Minimum depth to be passed to shapemapper2's "
                'make_reactivity_profile.py and render_figures.py.')
    parser.add_argument('--min-mapq', type=int, default=0,
            help="Minimum mapping quality to be passed to shapemapper2's "
                'shapemapper_mutation_parser in both of its --min_mapq and '
                '--min_qual arguments.')
    parser.add_argument('--max-bg', type=float, default=1,
            help='Maximum background mutation rate to be passed to '
                "shapemapper2's make_reactivity_profile.py and"
                'render_figures.py. Should be a float between 0 and 1.')
    parser.add_argument('--skip-plot', action='store_true', default=False,
            help='If specified, will not run render_figures.py (i.e. no '
                'reactivity plot will be produced.')
    parser.add_argument('--skip-shape', action='store_true', default=False,
            help='If specified, will not run tab_to_shape.py (i.e. no '
                'shape file will be produced.')

    parser.add_argument('--outdir', '-o', default='./',
            help='Directory to place output files within.')
    parser.add_argument('--keep', '-k', action='store_true', default=False,
            help='If specified, will not delete temp directory when done.')
