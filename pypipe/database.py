import sqlite3


class PipelineDatabase:

    def __init__(self, name):
        self.connection = sqlite3.connect(name)
        self.cursor = self.connection.cursor()

    def create_if_not_exists(self):
        self.cursor.execute('''create table if not exists files
        (file int, name text, format text, primary key (file))''')
        self.cursor.execute('''create table if not exists programs
        (program int, status int, name text, log text, out text, cmd text,
        primary key(program))''')
        self.cursor.execute('''create table if not exists files_programs
        (parent int, child int)''')
        self.cursor.execute('''create table if not exists programs_programs
        (parent int, child int, label text)''')
        self.cursor.execute('''create table if not exists programs_files
        (parent int, child int)''')

    def truncate_all(self):
        self.cursor.execute('delete from files')
        self.cursor.execute('delete from programs')
        self.cursor.execute('delete from files_programs')
        self.cursor.execute('delete from programs_programs')
        self.cursor.execute('delete from programs_files')

    def close(self):
        self.connection.close()

    def execute(self, sql):
        self.cursor.execute(sql)

    def fetchall(self):
        return self.cursor.fetchall()

    def commit(self):
        self.connection.commit()

