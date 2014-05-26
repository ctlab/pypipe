import os
import subprocess
import sys
import threading
import time
import tempfile

from pypipe import formats
from pypipe.paths import PYPIPE_DIR
import database
import abstractfile


COMPLETE, FAILED, RUNNING, NONE = range(4)


def add_arg_error_msg(name, type_, option=None):
    option = option and ('"' + option + '"') or 'argument'
    if issubclass(type_, abstractfile.File):
        t = type_.__name__.upper() + ' file'
    else:
        t = type_ + ' value'
    return '%s: %s must be %s' % (name, option, t)


def add_args_error_msg(name, type_, option=None):
    if option:
        msg = ': "%s" must be ' % option
    else:
        msg = ': expected list of '
    if issubclass(type_, abstractfile.File):
        t = type_.__name__.upper() + ' file'
    else:
        t = type_ + ' value'
    return name + msg + t


class PipelineNode:

    def __init__(self, name, log=None, out=None):
        self.name = name
        self.log = log or os.devnull
        self.out = out or self.log
        self.cmd = name.split(' ')
        self.children = set()
        self.parents = set()
        self.return_files = set()
        self.labels = {}
        self.count = 0
        self.number = None
        self.status = NONE
        self.bin_index = 0
        self.thread = None

    def join_pypipe_path(self):
        self.cmd[self.bin_index] = os.path.join(PYPIPE_DIR, self.cmd[self.bin_index])

    def thread_function(self):
        log = open(self.log, 'w')
        out = open(self.out, 'w')
        self.status = RUNNING
        try:
            status = subprocess.call(self.cmd, stdout=out, stderr=log)
        except OSError:
            self.join_pypipe_path()
            status = subprocess.call(self.cmd, stdout=out, stderr=log)
        out.close()
        log.close()
        msg = ' '.join(self.cmd) + (out.name != log.name and (
            ' > ' + out.name) or '')
        if status == 0:
            self.status = COMPLETE
            sys.stderr.write('"%s" complete\n' % msg)
        else:
            self.status = FAILED
            sys.stderr.write('"%s" failed\n' % msg)

    def run(self):
        self.thread = threading.Thread(target=self.thread_function, args=())
        self.thread.start()

    def reset(self):
        self.status = NONE
        for child in self.children:
            child.reset()

    def save(self, db):
        program = self.number
        status = self.status
        name = self.name
        log = self.log
        out = self.out
        cmd = self.cmd
        db.execute('''insert into programs (program, status, name, log, out, cmd) values
            (%d, %d, "%s", "%s", "%s", "%s")''' % (program, status, name, log, out, cmd))
        for c in self.return_files:
            child = c.number
            db.execute('insert into programs_files (parent, child) values (%d, %d)' % (program, child))
        for c in self.children:
            child = c.number
            label = c.labels[self]
            db.execute('''insert into programs_programs (parent, child, label)
                values (%d, %d, "%s")''' % (program, child, label))

    def add_arg(self, value, type_, option=None):
        if value is None or value is False:
            return
        if type(value) == int and type_ == float:
            value = float(value)
        if type(value) != type_:
            error_msg = add_arg_error_msg(self.name, type_, option)
            sys.exit(error_msg)
        if option:
            self.cmd.append(option)
        if isinstance(value, abstractfile.File):
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
            sys.exit(error_msg)
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
                if isinstance(v, abstractfile.File):
                    pipeline.add_edge(self, v)
            if issubclass(type_, abstractfile.File):
                value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)


class Pipeline:

    def __init__(self):
        self.to_run = set()
        self.running = []
        self.all_programs = []
        self.files = []

    def generate_to_run(self, node):
        if node.status == COMPLETE or node in self.to_run:
            return
        self.to_run.add(node)
        for child in node.children:
            child.count += 1
        for parent in node.parents:
            self.generate_to_run(parent)

    def remove_children_from_to_run(self, node):
        for child in node.children:
            try:
                self.to_run.remove(child)
            except KeyError:
                pass
            self.remove_children_from_to_run(child)

    def add_node(self, name, log=None, out=None, type_=None):
        program = PipelineNode(name, log, out)
        if type_ == 'jar':
            program.bin_index = 2
        self.all_programs.append(program)
        program.number = len(self.all_programs) - 1
        return program

    def add_file(self, f):
        if f not in self.files:
            self.files.append(f)
            f.number = len(self.files) - 1

    def add_edge(self, p, f):
        if f.program:
            p.parents.add(f.program)
            f.program.children.add(p)
            if f.program in p.labels:
                p.labels[f.program] += ', ' + f.get_name()
            else:
                p.labels[f.program] = f.get_name()
        else:
            f.next_programs.add(p)
            self.add_file(f)

    def run(self, i):
        self.generate_to_run(self.all_programs[i])
        while len(self.to_run) > 0 or len(self.running) > 0:
            for program in self.to_run:
                if program.count == 0:
                    self.running.append(program)
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
                    if program.status == FAILED:
                        self.remove_children_from_to_run(program)
            time.sleep(1)

    def reset(self, i):
        node = self.all_programs[i]
        node.reset()

    def reset_all(self):
        for node in self.all_programs:
            node.status = NONE

    def draw(self, img_name):
        graph = 'digraph {\n'
        i = 0
        for f in self.files:
            if f.next_programs:
                label = '%s\\n(%s)' % (f.get_name(), f.get_type())
                graph += '\t%d [label="%s" shape=rect];\n' % (i, label)
                f.graph_number = i
                i += 1
        for p in self.all_programs:
            graph += '\t%d [label="%s\\n(%d)" style=filled fillcolor=' % \
                     (i, p.name, p.number + 1)
            if p.status == COMPLETE:
                graph += '"green"'
            elif p.status == FAILED:
                graph += '"red"'
            elif p.status == RUNNING:
                graph += '"yellow"'
            else:
                graph += '"white"'
            graph += '];\n'
            p.graph_number = i
            i += 1
        for f in self.files:
            for p in f.next_programs:
                graph += '\t%d -> %d;\n' % (f.graph_number, p.graph_number)
        for p in self.all_programs:
            if len(p.children) > 0:
                for c in p.children:
                    label = c.labels[p]
                    graph += '\t%d -> %d [label="%s"];\n' % \
                             (p.graph_number, c.graph_number, label)
            else:
                for f in p.return_files:
                    label = '%s\\n(%s)' % (f.path, f.__class__.__name__)
                    graph += '\t%d [label="%s" shape=rect];\n' % (i, label)
                    graph += '\t%d -> %d;\n' % (p.graph_number, i)
                    i += 1
        graph += '}\n'
        dot_file = tempfile.mktemp()
        with open(dot_file, 'w+') as f:
            f.write(graph)
        try:
            subprocess.call(['dot', '-Tsvg', dot_file, '-o', img_name])
        except OSError:
            sys.stderr.write('Graphviz is not installed.\n')
        os.remove(dot_file)

    def save(self, db_name):
        db = database.PipelineDatabase(db_name)
        db.create_if_not_exists()
        db.truncate_all()
        for f in self.files:
            f.save(db)
        for p in self.all_programs:
            p.save(db)
        db.commit()
        db.close()

    def load(self, db_name):
        db = database.PipelineDatabase(db_name)
        db.execute('select * from programs')
        programs = []
        for program, status, name, log, out, cmd in db.fetchall():
            p = PipelineNode(name, log, out)
            p.number = program
            p.status = status
            p.cmd = eval(cmd)
            programs.append(p)
        self.all_programs = programs
        db.execute('select * from programs_programs')
        for parent, child, label in db.fetchall():
            p = self.all_programs[parent]
            c = self.all_programs[child]
            p.children.add(c)
            c.parents.add(p)
            c.labels[p] = label
        db.execute('select * from files')
        files = []
        for file_, name, format_ in db.fetchall():
            f = eval('formats.' + format_)(name, check=None)
            f.number = file_
            files.append(f)
        self.files = files
        db.execute('select * from files_programs')
        for parent, child in db.fetchall():
            p = self.files[parent]
            c = self.all_programs[child]
            p.next_programs.add(c)
        db.execute('select * from programs_files')
        for parent, child in db.fetchall():
            p = self.all_programs[parent]
            c = self.files[child]
            p.return_files.add(c)
        db.commit()
        db.close()

pipeline = Pipeline()

