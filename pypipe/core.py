import os
import subprocess
import sys
import threading
import time

from pypipe import formats
from pypipe.paths import PYPIPE_DIR, INSTALL_DIR


def add_arg_error_msg(name, type_, option=None):
    option = option and ('"' + option + '"') or 'argument'
    if issubclass(type_, formats._File):
        t = type_.__name__.upper() + ' file'
    else:
        t = type_ + ' value'
    return '%s: %s must be %s' % (name, option, t)


def add_args_error_msg(name, type_, options=None):
    if option:
        msg = ': "%s" must be ' % option
    else:
        msg = ': expected list of '
    if issubclass(type_, formats._File):
        t = type_.__name__.upper() + ' file'
    else:
        t = type_ + ' value'
    return name + msg + t


def join_pypipe_path(name, i):
    name = name.split(' ')
    name[i] = os.path.join(PYPIPE_DIR, name[i])
    name = ' '.join(name)
    return name


class Pipeline:

    def __init__(self):
        self.to_run = set()
        self.running = []
        self.all_programs = []
        self.input_files = []

    def generate_to_run(self, node):
        self.to_run.add(node)
        for parent in node.parents:
            self.generate_to_run(parent)

    def remove_children_from_to_run(self, node):
        for child in node.children:
            self.to_run.remove(child)
            self.remove_children_from_to_run(child)

    def add_node(self, name, log=None, out=None, type_=None):
        if type_ == 'jar':
            name = join_pypipe_path(name, 2)
        program = PipelineNode(name, log, out)
        self.all_programs.append(program)
        return program

    def add_edge(self, p, f):
        if f.program:
            p.parents.add(f.program)
            if p not in f.program.children:
                f.program.children.add(p)
                #f.program.labels...
                p.count += 1
            else:
                pass
                #labels...
        else:
            self.input_files.append(f)
            f.next_program = p

    def run(self, name, file_):
        try:
            self.generate_to_run(file_.program)
        except AttributeError:
            return
        while len(self.to_run) > 0 or len(self.running) > 0:
            for program in self.to_run:
                if program.count == 0:
                    program.run()
            for program in self.running:
                try:
                    self.to_run.remove(program)
                except KeyError:
                    pass
            running = [program for program in self.running]
            for program in running:
                if not program.thread.is_alive():
                    self.running.remove(program)
                    for child in program.children:
                        child.count -= 1
                    if program.ret != 0:
                        self.remove_children_from_to_run(program)
            time.sleep(1)


pipeline = Pipeline()


class PipelineNode:

    def __init__(self, name, log=None, out=None):
        if log and type(log) != str:
            sys.exit(name + ': "log" must be a string')
        if out and type(out) != str:
            sys.exit(name + ': "out" must be a string')
        self.name = name
        self.log = log or os.devnull
        self.out = out or self.log
        self.cmd = name.split(' ')
        self.children = set()
        self.parents = set()
        self.count = 0
        #self.labels = {}
        #self.return_files = []
        self.thread = threading.Thread(target=self.thread_function, args=())

    def thread_function(self):
        self.log = open(self.log, 'w')
        self.out = open(self.out, 'w')
        msg = ' '.join(self.cmd) + (self.out.name !=
                self.log.name and (' > ' + self.out.name) or '')
        try:
            self.ret = \
                subprocess.call(self.cmd, stdout=self.out, stderr=self.log)
        except OSError:
            self.cmd[0] = os.path.join(PYPIPE_DIR, self.cmd[0])
            self.ret = \
                subprocess.call(self.cmd, stdout=self.out, stderr=self.log)
        self.out.close()
        self.log.close()
        if self.ret == 0:
            sys.stderr.write('"%s" complete\n' % msg)
        else:
            sys.stderr.write('"%s" failed\n' % msg)

    def run(self):
        pipeline.running.append(self)
        self.thread.start()

    def add_arg(self, value, type_, option=None):
        if value is None or value == False:
            return
        if type(value) == int and type_ == float:
            value = float(value)
        if type(value) != type_:
            error_msg = add_arg_error_msg(self.name, type_, option)
            sys.exit(error_msg)
        if option:
            self.cmd.append(option)
        if isinstance(value, formats._File):
            self.cmd.append(value.path)
            pipeline.add_edge(self, value)
        elif type(value) != bool:
            self.cmd.append(str(value))

    def add_args(self, value, type_, min_len=1, delim=' ', option=None):
        if value is None:
            return
        if type(value) != list and type(value) != tuple:
            error_msg = add_args_error_msg(self.name, type_, option)
            sys.exit(error_msg)
        if len(value) < min_len:
            error_msg = '%s: list length must be >= %d' % (self.name, min_len)
        if not option:
            for v in value:
                self.add_arg(v, type_)
        else:
            for v in value:
                if type(v) == int and type_ == float:
                    v = float(v)
                if type(v) != type_:
                    error_msg = add_args_error_msg(self.name, type_, option)
                    sys.exit(error_msg)
                if isinstance(v, formats._File):
                    pipeline.add_edge(self, v)
            if issubclass(type_, formats._File):
                value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)

