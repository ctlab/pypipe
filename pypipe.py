#!/usr/bin/env python2

import argparse

from pypipe.formats import *
from pypipe.utils import generate_pipeline_graph, run_pipeline, reset_program


_parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
_parser.add_argument('pipeline', help='name of pipeline file', nargs=1,
        metavar="PIPELINE_FILE")
_group = _parser.add_mutually_exclusive_group(required=False)
_group.add_argument('--draw', action='store_true',
        help='draw pipeline to PNG')
_group.add_argument('--run', action='store', nargs=1,
        metavar='NODE_NAME', help='run pipeline')
_group.add_argument('--reset', action='store', nargs=1,
        metavar='NODE_NAME', help='reset pipeline from NODE')
_args = _parser.parse_args()

execfile(_args.pipeline)
if _args.draw:
    generate_pipeline_graph(_args.pipeline)
if _args.run:
    run_pipeline(_args.pipeline, eval(_args.run[0]))
if _args.reset:
    reset_program(_args.pipeline, eval(_args.reset[0]))

