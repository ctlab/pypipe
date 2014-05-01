#!/usr/bin/env python2

import argparse

from pypipe.formats import *
from pypipe.utils import run_pipeline, generate_pipeline_graph



_parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
_parser.add_argument('pipeline', help='name of pipeline file')
_parser.add_argument('--draw', action='store_true', 
        help='draw pipeline to PNG')
_parser.add_argument('--run', action='store', nargs=1, metavar='NODE',
        help='run pipeline')
_args = _parser.parse_args()

execfile(_args.pipeline)
if _args.draw:
    from pypipe.utils import generate_pipeline_graph
    generate_pipeline_graph(_args.pipeline)
if _args.run:
    from pypipe.utils import run_pipeline
    run_pipeline(_args.pipeline, eval(_args.run[0]))
