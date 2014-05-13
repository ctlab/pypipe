#!/usr/bin/env python2

import argparse

from pypipe.formats import *


_parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
_parser.add_argument('pipeline', help='name of pipeline file',
        metavar="PIPELINE_FILE")
_group = _parser.add_mutually_exclusive_group(required=True)
_group.add_argument('--draw', action='store_true',
        help='draw pipeline to PNG')
_group.add_argument('--run', action='store', nargs=1,
        metavar='NODE_NAME', help='run pipeline')
_group.add_argument('--reset', action='store', nargs=1,
        metavar='NODE_NAME', help='reset pipeline from NODE')
_group.add_argument('--resetall', action='store_true',
        help='reset pipeline from start noe')
_args = _parser.parse_args()

execfile(_args.pipeline)
if _args.draw:
    from pypipe.utils import generate_pipeline_graph
    generate_pipeline_graph(_args.pipeline)
elif _args.run:
    from pypipe.utils import run_pipeline
    run_pipeline(_args.pipeline, eval(_args.run[0]))
elif _args.reset:
    from pypipe.utils import reset_program
    reset_program(_args.pipeline, eval(_args.reset[0]))
elif _args.resetall:
    from pypipe.utils import reset_pipeline
    reset_pipeline(_args.pipeline)

