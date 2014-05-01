import os
import subprocess
import sys
import threading
import tempfile
import time

from pypipe import formats
from pypipe.paths import PYPIPE_DIR, INSTALL_DIR


_to_run = set()
_running = []
_all_programs = []
_input_files = []


def _add_arg_error_msg(name, type_, option=None):
    option = option and ("'" + option + "'") or "argument"
    msg = name + ": " + option + " must be "
    if issubclass(type_, formats._File):
        msg += type_.__name__.upper() + " file"
    elif type_ == int:
        msg += "an integer value"
    elif type_ == float:
        msg += "a float value"
    elif type == str:
        msg += "a string value"
    elif type_ == bool:
        msg += "a boolean value"
    return msg


def _add_args_error_msg(name, type_, options=None):
    msg = option and (": '%s' must be list of " % option) or \
            (": expected list of ")
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
    return msg


def _associate_program_with_file(p, f):
    if f.program:
        p.parents.add(f.program)
        if p not in f.program.children:
            f.program.children.add(p)
            f.program.labels[p] = f.path
            p.count += 1
        else:
            f.program.labels[p] += ", " + f.path
    else:
        _input_files.append(f)
        f.next_program = p


def _join_pypipe_path(name, i):
    name = name.split(" ")
    name[i] = os.path.join(PYPIPE_DIR, name[i])
    name = " ".join(name)
    return name


def _get_status_files(name):
    result = {}
    result["complete"] = name + ".complete"
    result["failed"] = name + ".failed"
    result["running"] = name + ".running"
    return result


def _read_status_file(filename):
    with open(filename, "r") as f:
        result = map(int, f.readline().split(" ")[:-1])
    return result


def _generate_to_run_list(status_files, node):
    _to_run.add(node)
    if len(node.parents) > 0:
        for parent in node.parents:
            _generate_to_run_list(status_files, parent)
    if os.path.exists(status_files["complete"]):
        complete = _read_status_file(status_files["complete"])
    else:
        complete = []
    for i in complete:
        try:
            _to_run.remove(_all_programs[i])
            for child in _all_programs[i].children:
                child.count -= 1
        except KeyError:
            pass


def _remove_children_from_to_run(node):
    try:
        _to_run.remove(node)
    except KeyError:
        pass
    if len(node.children) > 0:
        for child in node.children:
            _remove_children_from_to_run(child)


def _program_exists(program_name):
    paths = os.environ['PATH'].split(":")
    paths.append(PYPIPE_DIR)
    for path in paths:
        f = os.path.join(path, program_name)
        if os.path.isfile(f):
            return True
    return False


def _install_program(script_name, program_name):
    install_script = os.path.join(INSTALL_DIR, "install.sh")
    program_script = os.path.join(INSTALL_DIR, script_name)
    subprocess.call(["bash", install_script, program_script])


def _remove_children_numbers(node, program_to_index, numbers):
    try:
        numbers.remove(program_to_index[node])
    except ValueError:
        pass
    if len(node.children) > 0:
        for child in node.children:
            _remove_children_numbers(child, program_to_index, numbers)


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
        self.log = open(self.log, "w")
        self.output = open(self.output, "w")
        cmd_msg = " ".join(self.cmd) + (self.output.name !=
                self.log.name and (" > " + self.output.name) or "")
        try:
            self.ret = \
                subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        except OSError:
            self.cmd[0] = os.path.join(PYPIPE_DIR, self.cmd[0])
            self.ret = \
                subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        self.output.close()
        self.log.close()
        if self.ret == 0:
            sys.stderr.write("'%s' complete\n" % cmd_msg)
        else:
            sys.stderr.write("'%s' failed\n" % cmd_msg)

    def _run(self):
        _running.append(self)
        self.thread.start()

    def add_arg(self, value, type_, option=None):
        if value is None or value == False:
            return
        if type(value) == int and type_ == float:
            value = float(value)
        if type(value) != type_:
            error_msg = _add_arg_error_msg(self.name, type_, option)
            sys.exit(error_msg)
        if option:
            self.cmd.append(option)
        if isinstance(value, formats._File):
            self.cmd.append(value.path)
            _associate_program_with_file(self, value)
        elif type(value) != bool:
            self.cmd.append(str(value))

    def add_args(self, value, type_, min_len=1, delim=" ", option=None):
        if value is None:
            return
        if type(value) != list and type(value) != tuple:
            error_msg = _add_args_error_msg(self.name, type_, option)
            sys.exit(error_msg)
        if len(value) < min_len:
            error_msg = "%s: list length must be >= %d" % (self.name, min_len)
            sys.exit(error_msg)
        if not option:
            for v in value:
                self.add_arg(v, type_)
        else:
            for v in value:
                if type(v) == int and type_ == float:
                    v = float(v)
                if type(v) != type_:
                    error_msg = _add_args_error_msg(self.name, type_, option)
                    sys.exit(error_msg)
                if isinstance(v, formats._File):
                    _associate_program_with_file(self, v)
            if issubclass(type_, formats._File):
                value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)


def create_program(name, output=None, log=None, type_=None):
    if type_ == "jar":
        name = _join_pypipe_path(name, 2)
    program = _PipelineNode(name, output, log)
    _all_programs.append(program)
    return program


def install_program(script_name, program_name):
    if not _program_exists(program_name):
        print program_name, "is not installed. Installing..."
        _install_program(script_name)
        print program_name, "installed."


def run_pipeline(name, node):
    status_files = _get_status_files(name)
    i = 0
    program_to_index = {}
    for program in _all_programs:
        program_to_index[program] = i
        i += 1
    _generate_to_run_list(status_files, node.program)
    while len(_to_run) > 0 or len(_running) > 0:
        with open(status_files["running"], "w+") as f:
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
        running = [program for program in _running]
        for program in running:
            if not program.thread.is_alive():
                _running.remove(program)
                for child in program.children:
                    child.count -= 1
                if program.ret == 0:
                    log_file = status_files["complete"]
                else:
                    log_file = status_files["failed"]
                    _remove_children_from_to_run(program)
                with open(log_file, "a+") as f:
                    i = program_to_index[program]
                    f.write("%d " % i)
        time.sleep(1)
    os.remove(status_files["running"])


def generate_pipeline_graph(name):
    status_files = _get_status_files(name)
    dot_file = name + ".dot"
    png_file = name + ".png"
    complete = set()
    failed = set()
    running = set()
    if os.path.exists(status_files["complete"]):
        complete_i = _read_status_file(status_files["complete"])
        complete = set([_all_programs[i] for i in complete_i])
    if os.path.exists(status_files["failed"]):
        failed_i = _read_status_file(status_files["failed"])
        failed = set([_all_programs[i] for i in failed_i])
    if os.path.exists(status_files["running"]):
        running_i = _read_status_file(status_files["running"])
        running = set([_all_programs[i] for i in running_i])
    with open(dot_file, "w+") as f:
        d = {}
        i = 0
        f.write("digraph {\n")
        for e in _input_files:
            label = "%s\\n(%s)" % (e.path, e.__class__.__name__ )
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
            f.write('\t%d -> %d;\n' % (d[e], d[e.next_program]))
        for p in _all_programs:
            if len(p.children) > 0:
                for c in p.children:
                    f.write('\t%d -> %d [label="%s"];\n' % (d[p], d[c],
                        p.labels[c]))
            else:
                for out in p.return_files:
                    label = "%s\\n(%s)" % (out.path, out.__class__.__name__)
                    f.write('\t%d [label="%s" shape=rect];\n' % (i, label))
                    f.write('\t%d -> %d;\n' % (d[p], i))
                    i += 1
        f.write("}\n")
    try:
        subprocess.call(["dot", "-Tpng", dot_file, "-o", png_file])
    except OSError:
        sys.stderr.write("Graphviz is not installed\n")


def reset_program(name, node):
    status_files = _get_status_files(name)
    if os.path.exists(status_files["complete"]):
        complete = _read_status_file(status_files["complete"])
    program_to_index = {}
    i = 0
    for p in _all_programs:
        program_to_index[p] = i
        i += 1
    _remove_children_numbers(node.program, program_to_index, complete)
    with open(status_files["complete"], "w+") as f:
        for n in complete:
            f.write("%d " % n)

