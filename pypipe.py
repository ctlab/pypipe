#!/usr/bin/env python2

import argparse

from pypipe.formats import *


parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
parser.add_argument('pipeline', help='name of pipeline file',
        metavar="PIPELINE_FILE")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--draw', action='store',
        metavar='IMG_TYPE', help='draw pipeline to image file')
group.add_argument('--run', action='store', type=int,
        metavar='NODE_NUMBER', help='run pipeline')
group.add_argument('--reset', action='store', type=int,
        metavar='NODE_NUMBER', help='reset pipeline from node')
group.add_argument('--resetall', action='store_true',
        help='reset all pipeline nodes')
_args = parser.parse_args()

execfile(_args.pipeline)

from pypipe.core import pipeline
if _args.draw:
    pipeline.draw(_args.pipeline, _args.draw)
elif _args.run:
    pipeline.run(_args.pipeline, _args.run - 1)
elif _args.reset:
    pipeline.reset(_args.pipeline, _args.reset - 1)
elif _args.resetall:
    pipeline.reset_all(_args.pipeline)

