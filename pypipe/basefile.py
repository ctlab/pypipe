import ntpath
import os

import pypipe.baseexception


class FileNotExistsError(pypipe.baseexception.BaseException):

    def __init__(self, file_):
        if file_.suff:
            value = "One or more of files doesn't exist: " + ", ".join(file_.names())
        else:
            value = "File doesn't exist: " + file_.get_name()
        super(FileNotExistsError, self).__init__(value)


class File(object):

    def __init__(self, path, program, check=True, suff=None):
        self.program = program
        self.next_programs = set()
        self.path = path
        self.number = None
        self.suff = suff
        if program:
            self.program.return_files.add(self)
        elif check and not self.check():
            raise FileNotExistsError(self)

    def names(self):
        if self.suff:
            return [self.path + s for s in self.suff]
        return [self.path]

    def get_name(self):
        return ntpath.basename(self.path)

    def get_type(self):
        return self.__class__.__name__

    def check(self):
        for name in self.names():
            if not os.path.isfile(name):
                return False
        return True

    def save(self, db):
        file_ = self.number
        name = self.path
        format_ = self.get_type()
        suff = self.suff
        db.execute('insert into files (file, name, suff, format) values (%d, "%s", "%s", "%s")'
                   % (file_, name, suff, format_))
        for c in self.next_programs:
            child = c.number
            db.execute('insert into files_programs (parent, child) values (%d, %d)' % (file_, child))

