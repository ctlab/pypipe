#!/usr/bin/env python2

import argparse


parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
parser.add_argument('database', help='name of pipeline db',
        metavar="PIPELINE_DB_NAME")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--save', action='store',
        metavar='SCRIPT_NAME', help='save pipeline to db')
group.add_argument('--run', action='store', type=int,
        metavar='NODE_NUMBER', help='run pipeline')
group.add_argument('--reset', action='store', type=int,
        metavar='NODE_NUMBER', help='reset pipeline from node')
group.add_argument('--resetall', action='store_true',
        help='reset all pipeline nodes')
group.add_argument('--draw', action='store',
        metavar='IMG_NAME', help='draw pipeline to svg file')
group.add_argument('--renamefile', action='store', nargs=2,
        metavar=['1' '2'], help='draw pipeline to svg file')
_args = parser.parse_args()

from pypipe.core import pipeline
if _args.save:
    from pypipe.formats import *
    execfile(_args.save)
    pipeline.save(_args.database)
elif _args.run:
    pipeline.load(_args.database)
    pipeline.run(_args.run - 1)
    pipeline.save(_args.database)
elif _args.reset:
    pipeline.load(_args.database)
    pipeline.reset(_args.reset - 1)
    pipeline.save(_args.database)
elif _args.resetall:
    pipeline.load(_args.database)
    pipeline.reset_all()
    pipeline.save(_args.database)
elif _args.draw:
    pipeline.load(_args.database)
    pipeline.draw(_args.draw)
elif _args.renamefile:
    pipeline.load(_args.database)
    pipeline.rename_file()
    pipeline.save(_args.database)

