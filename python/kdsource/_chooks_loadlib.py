_loadlib_cache = [None,None]
def _loadlib_kdsource():
    if _loadlib_cache[1]:
        return _loadlib_cache[1]
    from .config import config
    _loadlib_cache[0] = _loadlib_frompath(_get_mcpl_shlibpath())
    _loadlib_cache[1] = _loadlib_frompath(config('shlibpath'))
    return _loadlib_cache[1]

def _loadlib_mcpl():
    if _loadlib_cache[0] is None:
        _loadlib_kdsource()
    return _loadlib_cache[0]

def _get_mcpl_shlibpath():
    import subprocess
    rv = subprocess.run( ['mcpl-config','--show','shlibpath'],
                         check = True, capture_output = True )
    if rv.returncode or rv.stderr:
        raise RuntimeError('Problems invoking mcpl-config for shlibpath')
    import pathlib
    return pathlib.Path(rv.stdout.decode().strip()).absolute().resolve()

def _loadlib_frompath(libpath):
    import ctypes
    try:
        lib = ctypes.CDLL(libpath)
    except TypeError:
        lib = None
    if lib is None:
        #On Windows, providing a Path rather than str can yield TypeError:
        lib = ctypes.CDLL(str(libpath))
    return lib

