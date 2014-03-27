import os
import subprocess
import sys
import threading
import time

import formats


_to_run = []
_running = []


class _PipelineNode:

    def __init__(self, name, output=None):
        if output and type(output) != str:
            sys.exit(name + ": output must be string")
        self.name = name
        self.cmd = name.split(" ")
        self.children = set()
        self.count = 0
        self.output = output
        self.thread = threading.Thread(target=self.thread_function, args=())
        
    def thread_function(self):
        msg = " ".join(self.cmd)+(self.output and (" > " + self.output) or "")
        sys.stderr.write(msg + "\n")
        FNULL = open(os.devnull, "w")
        if self.output:
            with open(self.output, "w") as output_file:
                subprocess.call(self.cmd, stdout=output_file, stderr=FNULL)
        else:
            subprocess.call(self.cmd, stdout=FNULL, stderr=FNULL)
        FNULL.close()

    def add_arg(self, value, type_, option=None):
        if not value:
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
        if not value:
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
                    value = [v.path for v in value]
            self.add_arg(delim.join(map(str, value)), str, option)

    def run(self):
        _running.append(self)
        self.thread.start()


def create_program(name, output=None):
    program = _PipelineNode(name, output)
    _to_run.append(program)
    return program


def run_pipeline():
    while len(_to_run) > 0:
        for program in _to_run:
            if program.count == 0:
                program.run()
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

