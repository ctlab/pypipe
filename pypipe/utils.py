import os
import subprocess
import sys
import threading
import time
import tempfile

from pypipe import formats


_pypipe_dir = os.path.join(os.environ['HOME'], '.pypipe')
_install_dir = os.path.join(_pypipe_dir, "install-scripts")
_complete_file = "example2.py.complete"
_failed_file = "example2.py.failed"
_running_file = "example2.py.running"

_to_run = set()
_running = []
_all_programs = []
_input_files = []


class _PipelineNode:

    def __init__(self, name, log=None, output=None):
        if output and type(output) != str:
            sys.exit(name + ": output must be string")
        if log and type(log) != str:
            sys.exit(name + ": log must be string")
        self.name = name
        self.cmd = name.split(" ")
        self.children = set()
        self.parents = set()
        self.count = 0
        self.labels = {}
        self.return_files = []
        if log:
            self.log = log
        else:
            self.log = os.devnull
        if output:
            self.output = output
        else:
            self.output = self.log
        self.thread = threading.Thread(target=self._thread_function, args=())
        
    def _thread_function(self):
        cmd_msg = " ".join(self.cmd) + (self.output.name !=
                self.log.name and (" > " + self.output.name) or "")
        #sys.stderr.write(cmd_msg + "\n")  # debug
        try:
            self.ret = \
                subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        except OSError:
            self.cmd[0] = os.path.join(_pypipe_dir, self.cmd[0])
            self.ret = \
                subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        if self.ret != 0:
            for child in self.children:
                _to_run.remove(child)
        self.output.close()
        self.log.close()
        if self.ret == 0:
            sys.stderr.write("'" + cmd_msg + "' complete\n")  # debug
        else:
            sys.stderr.write("'" + cmd_msg + "' failed\n")  # debug

    def add_arg(self, value, type_, option=None):
        if value is None or value == False:
            return
        if type(value) == int and type_ == float:
            value = float(value)

        if option:
            option2 = "'" + option + "'"
        else:
            option2 = "argument"
        msg = self.name + ": " + option2 + " must be "
        if issubclass(type_, formats._File):
            msg += type_.__name__.upper() + " file"
        elif type_ == int:
            msg += "an integer value"
        elif type_ == float:
            msg += "a float value"
        elif type_ == str:
            msg += "a string value"
        elif type_ == bool:
            msg += "a boolean value"

        if type(value) != type_:
            sys.exit(msg)
        if option:
            self.cmd.append(option)
        if isinstance(value, formats._File):
            self.cmd.append(value.path)
            if value.program: 
                self.parents.add(value.program)
                if self not in value.program.children:
                    value.program.children.add(self)
                    value.program.labels[self] = value.path
                    self.count += 1
                else:
                    value.program.labels[self] += ", " + value.path
            else:
                _input_files.append(value)
                value.next_program = self
        elif type(value) != bool:
            self.cmd.append(str(value))

    def add_args(self, value, type_, min_len=1, delim=" ", option=None):
        if value is None:
            return

        if option:
            msg = self.name + ": '" + option + "' must be list of "
        else:
            msg = self.name + ": expected list of "
        if issubclass(type_, formats._File):
            msg += type_.__name__.upper() + " files"
        elif type_ == int:
            msg += "integer values"
        elif type_ == float:
            msg += "float values"
        elif type_ == str:
            msg += "string values"
        elif type_ == bool:
            msg += "boolean values"

        if type(value) != list and type(value) != tuple:
            sys.exit(msg)
        if len(value) < min_len:
            sys.exit(self.name + ": list length must be >= " + str(min_len))

        if not option:
            for v in value:
                self.add_arg(v, type_)
        else:
            for v in value:
                if type(v) == int and type_ == float:
                    v = float(v)
                if type(v) != type_:
                    sys.exit(msg)
                if isinstance(v, formats._File):
                    if v.program:
                        self.parents.add(v.program)
                        if self not in v.program.children:
                            v.program.children.add(self)
                            v.program.labels[self] = v.path
                            self.count += 1
                        else:
                            v.program.labels[self] += ", " + v.path
                    else:
                        _input_files.append(v)
                        v.next_program = self
            if issubclass(type_, formats._File):
                value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)

    def _run(self):
        self.log = open(self.log, "w")
        self.output = open(self.output, "w")
        _running.append(self)
        self.thread.start()


def _join_pypipe_path(name, i):
    name = name.split(" ")
    name[i] = os.path.join(_pypipe_dir, name[i])
    name = " ".join(name)
    return name


def create_program(name, output=None, log=None, type_=None):
    if type_ == "jar":
        name = _join_pypipe_path(name, 2)
    program = _PipelineNode(name, output, log)
    _all_programs.append(program)
    return program


def _append_to_run(node):
    _to_run.add(node)
    if len(node.parents) > 0:
        for parent in node.parents:
            _generate_to_run(parent)
    if os.path.exists(_complete_file):
        with open(_complete_file, "r") as f:
            complete = map(int, f.readline().split(" ")[:-1])
        for i in complete:
            try:
                _to_run.remove(_all_programs[i])
                for child in _all_programs[i].children:
                    child.count -= 1
            except KeyError:
                pass


def _remove_to_run(node):
    try:
        _to_run.remove(node)
    except KeyError:
        pass
    if len(node.children) > 0:
        for child in node.children:
            _remove_to_run(child)


def run_pipeline(node):
    program_to_index = {}
    i = 0
    for program in _all_programs:
        program_to_index[program] = i
        i += 1
    _append_to_run(node.program)
    while len(_to_run) > 0 or len(_running) > 0:
        with open(_running_file, "w+") as f:
            for program in _running:
                f.write("%d " % program_to_index[program])
        for program in _to_run:
            if program.count == 0:
                program._run()
        for program in _running:
            try:
                _to_run.remove(program)
            except KeyError:
                pass
        running = [p for p in _running]
        for program in running:
            if not program.thread.is_alive():
                if program.ret == 0:
                    log_file = _complete_file
                else:
                    log_file = _failed_file
                _running.remove(program)
                program.log.close()
                program.output.close()
                with open(log_file, "a+") as f:
                    i = program_to_index[program]
                    f.write("%d " % i)
                for child in program.children:
                    child.count -= 1
        time.sleep(1)
    os.remove(_running_file)


def _program_exists(program_name):
    paths = os.environ['PATH'].split(":")
    paths.append(_pypipe_dir)
    for path in paths:
        f = os.path.join(path, program_name)
        if os.path.isfile(f):
            return True
    return False


def _install_program(script_name):
    install_script = os.path.join(_install_dir, "install.sh")
    program_script = os.path.join(_install_dir, script_name)
    subprocess.call(["bash", install_script, program_script])


def install_program(script_name, program_name):
    if not _program_exists(program_name):
        print program_name, "is not installed. Installing..."
        _install_program(script_name)
        print program_name, "installed."


def generate_pipeline_graph(filename):
    dot_name = filename + ".dot"
    png_name = filename + ".png"
    complete = set()
    failed = set()
    running = set()
    if os.path.exists(_complete_file):
        with open(_complete_file, "r") as f:
            complete_i = map(int, f.readline().split(" ")[:-1])
            complete = set([_all_programs[i] for i in complete_i])
    if os.path.exists(_failed_file):
        with open(_failed_file, "r") as f:
            failed_i = map(int, f.readline().split(" ")[:-1])
            failed = set([_all_programs[i] for i in failed_i])
    if os.path.exists(_running_file):
        with open(_running_file, "r") as f:
            running_i = map(int, f.readline().split(" ")[:-1])
            running = set([_all_programs[i] for i in running_i])
    with open(dot_name, "w+") as f:
        d = {}
        i = 0
        f.write("digraph {\n")
        for e in _input_files:
            label = e.path + "\\n(" + e.__class__.__name__ + ")"
            f.write('\t%d [label="%s" shape=rect];\n' % (i, label))
            d[e] = i
            i += 1
        for p in _all_programs:
            f.write('\t%d [label="%s"' % (i, p.name))
            if p in complete:
                f.write(' style=filled fillcolor="green"')
            elif p in failed:
                f.write(' style=filled fillcolor="red"')
            elif p in running:
                f.write(' style=filled fillcolor="yellow"')
            f.write(']\n')
            d[p] = i
            i += 1
        for e in _input_files:
            f.write('\t%d -> %d;\n' % (d[e] , d[e.next_program]))
        for p in _all_programs:
            if len(p.children) > 0:
                for c in p.children:
                    f.write('\t%d -> %d [label="%s"];\n' % (d[p] , d[c], 
                        p.labels[c]))
            else:
                for out in p.return_files:
                    label = out.path + "\\n(" + out.__class__.__name__ + ")"
                    f.write('\t%d [label="%s" shape=rect];\n' % (i, label))
                    f.write('\t%d -> %d;\n' % (d[p] , i))
                    i += 1
        f.write("}\n")
    try:
        subprocess.call(["dot", "-Tpng", dot_name, "-o", png_name])
    except OSError:
        sys.stderr.write("Graphviz is not installed\n")

