import os
import subprocess
import sys
import threading
import time
import ntpath
import sqlite3

from pypipe import formats
from pypipe.paths import PYPIPE_DIR, INSTALL_DIR


RUNNING, COMPLETE, FAILED, NONE = range(4)


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


class Database:

    def __init__(self, name):
        self.connection = sqlite3.connect(name + '.db')
        self.cursor = self.connection.cursor()
        self.cursor.execute('''create table if not exists nodes
                (node int, status int, primary key (node))''')
        self.connection.commit()

    def close(self):
        self.connection.close()

    def save_node_status(self, node, status):
        self.cursor.execute('''insert or replace into nodes (node, status)
                   values (%d, %d)''' % (node.number, status))
        self.connection.commit()

    def _remove_node_statuses(self, node):
        self.cursor.execute('delete from nodes where node = %d' % node.number)
        for child in node.children:
            self._remove_node_statuses(child)

    def remove_node_statuses(self, node):
        self._remove_node_statuses(node)
        self.connection.commit()

    def remove_all(self):
        self.cursor.execute('delete from nodes')
        self.connection.commit()

    def remove_running_nodes(self):
        self.cursor.execute('''update nodes set status = %d
                where status = %d''' % (FAILED, RUNNING))

    def get_nodes_status(self):
        self.cursor.execute('select * from nodes')
        return self.cursor.fetchall()


class Pipeline:

    def __init__(self):
        self.to_run = set()
        self.running = []
        self.all_programs = []
        self.input_files = []

    def generate_to_run(self, node, statuses):
        n = node.number
        if n in statuses and statuses[n] == COMPLETE:
            return
        self.to_run.add(node)
        for child in node.children:
            child.count += 1
        for parent in node.parents:
            self.generate_to_run(parent, statuses)

    def remove_children_from_to_run(self, node):
        for child in node.children:
            self.to_run.remove(child)
            self.remove_children_from_to_run(child)

    def add_node(self, name, log=None, out=None, type_=None):
        if type_ == 'jar':
            name = join_pypipe_path(name, 2)
        program = PipelineNode(name, log, out)
        self.all_programs.append(program)
        program.number = len(self.all_programs) - 1
        return program

    def add_edge(self, p, f):
        if f.program:
            p.parents.add(f.program)
            f.program.children.add(p)
            if f.program in p.input_files:
                p.input_files[f.program].append(f)
            else:
                p.input_files[f.program] = [f]
        else:
            self.input_files.append(f)
            f.next_programs.append(p)

    def reset(self, name, i):
        db = Database(name)
        db.remove_node_statuses(self.all_programs[i])
        db.close()

    def reset_all(self, name):
        db = Database(name)
        db.remove_all()
        db.close()

    def draw(self, name, img_type):
        db = Database(name)
        statuses = dict(db.get_nodes_status())
        db.close()
        graph = 'digraph {\n'
        i = 0
        for f in self.input_files:
            label = '%s\\n(%s)' % \
                (ntpath.basename(f.path), f.__class__.__name__)
            graph += '\t%d [label="%s" shape=rect];\n' % (i, label)
            f.number = i
            i += 1
        for p in self.all_programs:
            graph += '\t%d [label="%s\\n(%d)" style=filled fillcolor=' % \
                    (i, p.name, p.number + 1)
            if p.number in statuses:
                status = statuses[p.number]
            else:
                status = NONE
            if status == COMPLETE:
                graph += '"green"'
            elif status == FAILED:
                graph += '"red"'
            elif status == RUNNING:
                graph += '"yellow"'
            else:
                graph += '"white"'
            graph += '];\n'
            p.graph_number = i
            i += 1
        for f in self.input_files:
            for p in f.next_programs:
                graph += '\t%d -> %d;\n' % (f.number, p.graph_number)
        for p in self.all_programs:
            if len(p.children) > 0:
                for c in p.children:
                    label = ','.join([f.path for f in c.input_files[p]])
                    graph += '\t%d -> %d [label="%s"];\n' % \
                            (p.graph_number, c.graph_number, label)
            else:
                for f in p.return_files:
                    label = '%s\\n(%s)' % (f.path, f.__class__.__name__)
                    graph += '\t%d [label="%s" shape=rect];\n' % (i, label)
                    graph += '\t%d -> %d;\n' % (p.graph_number, i)
                    i += 1
        graph += '}\n'
        dot_file = name + '.dot'
        img_file = name + '.' + img_type
        with open(dot_file, 'w+') as f:
            f.write(graph)
        try:
            subprocess.call(['dot', '-T' + img_type, dot_file, '-o', img_file])
        except OSError:
            sys.stderr.write('Graphviz is not installed.\n')


    def run(self, name, i):
        db = Database(name)
        statuses = dict(db.get_nodes_status())
        self.generate_to_run(self.all_programs[i], statuses)
        while len(self.to_run) > 0 or len(self.running) > 0:
            for program in self.to_run:
                if program.count == 0:
                    db.save_node_status(program, RUNNING)
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
                        db.save_node_status(program, FAILED)
                    else:
                        db.save_node_status(program, COMPLETE)
            time.sleep(1)
        db.remove_running_nodes()
        db.close()


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
        self.input_files = {}
        self.return_files = []
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

