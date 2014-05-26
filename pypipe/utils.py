import subprocess

from pypipe.core import pipeline 
from pypipe.paths import PYPIPE_DIR, INSTALL_SCRIPTS_DIR
from pypipe.formats import *


def _program_exists(program_name):
    paths = os.environ['PATH'].split(':')
    paths.append(PYPIPE_DIR)
    for path in paths:
        f = os.path.join(path, program_name)
        if os.path.isfile(f):
            return True
    return False


def install_program(script_name, program_name):
    if not _program_exists(program_name):
        print program_name, 'is not installed. Installing...'
        install_scripts_path = os.path.join(os.getcwd(), INSTALL_SCRIPTS_DIR)
        install_script = os.path.join(install_scripts_path, 'install.sh')
        program_script = os.path.join(install_scripts_path, script_name)
        subprocess.call(['sh', install_script, program_script])
        print program_name, 'installed.'


def _handle_type(t, kwargs):
    if type(t) == dict:
        ok = False
        for k in t:
            kk = k.replace('-', '_')
            if kk in kwargs and kwargs[kk]:
                ok = True
                type_ = t[k]
        if not ok:
            type_ = t['']
    else:
        type_ = t
    return type_


def _handle_arg(cmd, o, t, kwargs, option_none=False):
    option = o.replace('*', '')
    key = option.replace('-', '_')
    type_ = _handle_type(t, kwargs)
    mandatory = (o[-1] == '*') and True or False
    value = (key in kwargs) and kwargs[key] or None
    if mandatory and value is None:
        sys.exit('%s: %s is a mandatory argument' % (cmd, key))
    if option_none:
        option = None
    return {
        'option': option,
        'type': type_,
        'key': key,
        'value': value,
    }


def _handle_return(ret, kwargs):
    result = []
    for r in ret:
        key = r['arg'].replace('-', '_')
        name = kwargs[key] + r['suffix']
        type_ = _handle_type(r['type'], kwargs)
        if type_ != None:
            result.append({
                'key': key,
                'name': name,
                'type': type_,
            })
    return result


def tool(config):
    def wrapper(**kwargs):
        conf = config()
        cmd = conf['cmd']
        type_ = conf['type']

        log_k = conf['log']
        log = log_k in kwargs and kwargs[log_k] or None
        if log and type(log) != str:
            sys.exit('%s: log argument "%s" must be a string' % (cmd, log_k))

        args = []
        named_args = conf['args']['named']
        for k in named_args:
            arg = _handle_arg(cmd, k, named_args[k], kwargs)
            args.append(arg)
        unnamed_args = conf['args']['unnamed']
        for o, t in unnamed_args:
            arg = _handle_arg(cmd, o, t, kwargs, option_none=True)
            args.append(arg)

        return_info = _handle_return(conf['out']['return'], kwargs)
        out_keys = []
        if conf['out']['redirect']:
            out = return_info[0]['name']
            for r in return_info:
                out_keys.append(r['key'])
        else:
            out = None


        allowable_keys = set([a['key'] for a in args])
        allowable_keys.add(log_k)
        for k in kwargs:
            if k not in allowable_keys:
                sys.exit('%s: unexpected key "%s"' % (cmd, k))
        
        program = pipeline.add_node(cmd, log, out, type_)
        for arg in args:
            if arg['key'] in out_keys:
                continue
            if type(arg['type']) == list:
                type_ = arg['type'][0]
                min_len = arg['type'][1]
                delim = arg['type'][2]
                value = arg['value']
                option = arg['option']
                program.add_args(value, type_, min_len, delim, option)
            else:
                program.add_arg(arg['value'], arg['type'], arg['option'])

        return_values = []
        for r in return_info:
            name = r['name']
            init = r['type']
            file_ = init(name, program)
            return_values.append(file_)
            pipeline.add_file(file_)
        if len(return_values) == 1:
            return return_values[0]
        return return_values

    return wrapper

