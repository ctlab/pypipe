import os
import subprocess
import sys
import threading
import time
import tempfile

from pypipe import formats

_pypipe_dir = os.path.join(os.environ['HOME'], '.pypipe')
_install_dir = os.path.join(_pypipe_dir, "install-scripts")

_to_run = []
_running = []


class _PipelineNode:

    def __init__(self, name, log=None, output=None):
        if output and type(output) != str:
            sys.exit(name + ": output must be string")
        if log and type(log) != str:
            sys.exit(name + ": log must be string")
        self.name = name
        self.cmd = name.split(" ")
        self.children = set()
        self.count = 0
        if log:
            self.log = open(log, "w")
        else:
            self.log = open(os.devnull, "w")
        if output:
            self.output = open(output, "w")
        else:
            self.output = self.log
        self.thread = threading.Thread(target=self._thread_function, args=())
        
    def _thread_function(self):
        cmd_msg = " ".join(self.cmd) + (self.output.name !=
                self.log.name and (" > " + self.output.name) or "")
        sys.stderr.write(cmd_msg + "\n")  # debug
        try:
            subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        except OSError:
            self.cmd[0] = os.path.join(_pypipe_dir, self.cmd[0])
            subprocess.call(self.cmd, stdout=self.output, stderr=self.log)
        self.output.close()
        self.log.close()
        sys.stderr.write("'" + cmd_msg + "' complete\n")  # debug

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
            if value.program and self not in value.program.children:
                value.program.children.add(self)
                self.count += 1
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
                    if v.program and self not in v.program.children:
                        v.program.children.add(self)
                        self.count += 1
            if issubclass(type_, formats._File):
                value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)

    def _run(self):
        _running.append(self)
        self.thread.start()


def create_program(name, output=None, log=None):
    program = _PipelineNode(name, output, log)
    _to_run.append(program)
    return program


def run_pipeline():
    while len(_to_run) > 0:
        for program in _to_run:
            if program.count == 0:
                program._run()
        for program in _running:
            try:
                _to_run.remove(program)
            except ValueError:
                pass
        for program in _running:
            if not program.thread.is_alive():
                _running.remove(program)
                for child in program.children:
                    child.count -= 1
        time.sleep(1)


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

