import sys
import ntpath
import os


class File(object):

    def __init__(self, path, program, check=True, suff=[]):
        self.program = program
        self.next_programs = set()
        self.path = path
        self.number = None
        if suff:
            self.names = [path + s for s in suff]
        else:
            self.names = [path]
        if program:
            self.program.return_files.add(self)
        elif check and not self.check():
            if suff:
                sys.exit("One or more of files doesn't exist: " + ", ".join(self.names))
            else:
                sys.exit("File doesn't exist: " + self.get_name())

    def get_name(self):
        return ntpath.basename(self.path)

    def get_type(self):
        return self.__class__.__name__

    def check(self):
        for name in self.names:
            if not os.path.isfile(name):
                return False
        return True

    def save(self, db):
        file_ = self.number
        name = self.path
        format_ = self.get_type()
        db.execute('insert into files (file, name, format) values (%d, "%s", "%s")' % (file_, name, format_))
        for c in self.next_programs:
            child = c.number
            db.execute('insert into files_programs (parent, child) values (%d, %d)' % (file_, child))

