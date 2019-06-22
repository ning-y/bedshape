#! /usr/bin/env python3

import os

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

MUT_PARSER_BIN = os.path.join(THIS_DIR, 'shapemapper_mutation_parser')
MUT_COUNTER_BIN = os.path.join(THIS_DIR, 'shapemapper_mutation_counter')

PROFILER_BIN = os.path.join(THIS_DIR, 'make_reactivity_profiles.py')
NORMALISER_BIN = os.path.join(THIS_DIR, 'normalize_profiles.py')
RENDERER_BIN = os.path.join(THIS_DIR, 'render_figures.py')
TAB2SHAPE_BIN = os.path.join(THIS_DIR, 'tab_to_shape.py')
