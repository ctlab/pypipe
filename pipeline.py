import subprocess
import threading
import time


to_run = []
running = []


class PipelineNode:

    def __init__(self, cmd, output=None):
        self.cmd = cmd
        self.children = set()
        self.count = 0
        self.output = output
        self.thread = threading.Thread(target=self.thread_function, args=())
        
    def thread_function(self):
        if self.output:
            with open(self.output, "w") as output_file:
                subprocess.call(self.cmd, stdout=output_file)
        else:
            subprocess.call(self.cmd)

    def run(self):
        running.append(self)
        self.thread.start()


def create_program(cmd, files, output=None):
    program = PipelineNode(cmd, output)
    to_run.append(program)
    for f in files:
        if f.program and program not in f.program.children:
            f.program.children.add(program)
            program.count += 1
    return program


def run_pipeline():
    while len(to_run) > 0:
        for program in to_run:
            if program.count == 0:
                program.run()
        for program in running:
            try:
                to_run.remove(program)
            except ValueError:
                pass
        for program in running:
            if not program.thread.is_alive():
                running.remove(program)
                for child in program.children:
                    child.count -= 1
        time.sleep(1)

