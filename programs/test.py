import formats
from pipeline import create_program


def touch(output):
    cmd = ["touch", output]
    program = create_program(cmd, [])
    return formats.Sam(output, program)


def cp(f, g):
    cmd = ["cp", f.path, f.path + "-copy"]
    program = create_program(cmd, [f, g])
    return formats.Sam(f.path + "-copy", program)


def sleep():
    cmd = ["sleep", "20"]
    program = create_program(cmd, [])
    return formats.Sam("20", program)
