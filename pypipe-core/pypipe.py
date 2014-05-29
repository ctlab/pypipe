#!/usr/bin/env python2

import argparse


parser = argparse.ArgumentParser(
        description="Bioinformatics pipelines framework")
subparsers = parser.add_subparsers(dest='subparser')

create = subparsers.add_parser('create',
        help='create pipeline file')
create.add_argument('script', metavar='SCRIPT',
        help='pypipe script file')
create.add_argument('pipeline', metavar='PIPELINE',
        help='pipeline file that will be created')

run = subparsers.add_parser('run',
        help='run pipeline')
run.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')
run.add_argument('n', metavar='N', type=int,
        help='program node number that you need to run')

run_all = subparsers.add_parser('run-all',
        help='run all nodes in pipeline')
run_all.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')

reset = subparsers.add_parser('reset',
        help='reset pipeline node')
reset.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')
reset.add_argument('n', metavar='N', type=int,
        help='program node number that you need to reset')

reset_all = subparsers.add_parser('reset-all',
        help='reset all pipeline nodes')
reset_all.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')

draw = subparsers.add_parser('draw',
        help='draw pipeline to svg image file')
draw.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')
draw.add_argument('img', metavar='IMG',
        help='name of image that will be created')

rename = subparsers.add_parser('rename-file',
        help='rename input file')
rename.add_argument('pipeline', metavar="PIPELINE",
        help='pipeline file')
rename.add_argument('n', metavar='N', type=int,
        help='file node number that you need to rename')
rename.add_argument('new_name', metavar='NEW_NAME', 
        help='new name of file')

_args = parser.parse_args()

from pypipe.core import pipeline
if _args.subparser == 'create':
    from pypipe.formats import *
    try:
        execfile(_args.script)
    except Exception as e:
        print e
    pipeline.save(_args.pipeline)
elif _args.subparser == 'run':
    try:
        pipeline.load(_args.pipeline)
        pipeline.run(_args.n - 1)
    except Exception as e:
        print e
    pipeline.save(_args.pipeline)
elif _args.subparser == 'run-all':
    try:
        pipeline.load(_args.pipeline)
        pipeline.run()
    except Exception as e:
        print e
    pipeline.save(_args.pipeline)
elif _args.subparser == 'reset':
    try:
        pipeline.load(_args.pipeline)
        pipeline.reset(_args.n - 1)
    except Exception as e:
        print e
    pipeline.save(_args.pipeline)
elif _args.subparser == 'reset-all':
    try:
        pipeline.load(_args.pipeline)
    except Exception as e:
        print e
    pipeline.reset_all()
    pipeline.save(_args.pipeline)
elif _args.subparser == 'draw':
    try:
        pipeline.load(_args.pipeline)
        pipeline.draw(_args.img)
    except Exception as e:
        print e
elif _args.subparser == 'change-file':
    try:
        pipeline.load(_args.pipeline)
        pipeline.change_file(_args.n - 1, _args.new_name)
    except Exception as e:
        print e
    pipeline.save(_args.pipeline)

