def _str2cstr(s):
    #converts any string (str,bytes,unicode,path) to bytes
    if hasattr(s,'__fspath__'):
        s = str(s)
    try:
        return s if isinstance(s,bytes) else s.encode('utf8')
    except UnicodeEncodeError as e:
        from .exceptions import NCBadInput
        raise NCBadInput("Only unicode strings are supported") from e

def _cstr2str(s):
    #converts bytes object to str
    try:
        return s if isinstance(s,str) else s.decode('utf8')
    except UnicodeDecodeError as e:
        from .exceptions import NCBadInput
        raise NCBadInput("Only UTF8-encoded C-strings are supported") from e

__cache_fcts=[None]
def _load_fcts():
    if __cache_fcts[0] is not None:
        return __cache_fcts[0]
    from ._chooks_loadlib import _loadlib_kdsource
    import ctypes
    allfcts = {}
    lib = _loadlib_kdsource()

    ###### resample_to_mcpl
    kds_rng_fct_t = ctypes.CFUNCTYPE( ctypes.c_double )
    rawfct_rs = lib.kdsource_resample_to_mcpl
    rawfct_rs.restype=None
    rawfct_rs.argtypes=( kds_rng_fct_t,
                         ctypes.c_char_p,
                         ctypes.c_char_p,
                         ctypes.c_uint64 )
    def fct_rs( rng_fct, kds_sourcefile, destination_mcpl, nout ):
        rawfct_rs( kds_rng_fct_t(rng_fct),
                   _str2cstr(kds_sourcefile),
                   _str2cstr(destination_mcpl),
                   ctypes.c_uint64(nout) )
    allfcts['resample_to_mcpl'] = fct_rs

    ###### write_mcpl
    dblptr = ctypes.POINTER(ctypes.c_double)
    int32ptr = ctypes.POINTER(ctypes.c_int32)
    uint32ptr = ctypes.POINTER(ctypes.c_uint32)
    rawfct_wm = lib.kdsource_write_mcpl
    rawfct_wm.restype = None
    rawfct_wm.argtypes=( ctypes.c_char_p,ctypes.c_uint64,ctypes.c_uint,
                         dblptr,dblptr,dblptr,
                         dblptr,dblptr,dblptr,dblptr,dblptr,dblptr,int32ptr,
                         dblptr,dblptr,dblptr, uint32ptr)
    def get_ndarray_to_cptr():
        import numpy as np
        def ndarray_fix_type(arr, dtypename):
            a = np.asarray(arr)
            dt = np.dtype(dtypename)
            if a.dtype == dt and a.flags['C_CONTIGUOUS']:
                return a
            return a.astype(dt, copy=True)
        def ndarray_to_cptr(keepalive,a,dtypename='float64',cptr=dblptr):
            keepalive.append( ndarray_fix_type(a,dtypename) )
            return keepalive[-1].ctypes.data_as(cptr)
        return ndarray_to_cptr
    def fct_wm( filename, nparticles,
                ekin, x, y, z, ux, uy, uz, time, weight, pdgcode,
                polx, poly, polz, userflags, double_prec ):
        import numbers
        flags = 1 if double_prec else 0
        keepalive=[]
        a2c = get_ndarray_to_cptr()
        if isinstance(pdgcode, numbers.Integral):
            flags += 2
            pdgcode_array =  (ctypes.c_int32 * 1)(int(pdgcode))
            pdgcode_arg = ctypes.cast(pdgcode_array, int32ptr)
        else:
            pdgcode_arg = a2c(keepalive,pdgcode,'int32',int32ptr)
        if isinstance(weight, numbers.Number):
            flags += 4
            weight_array =  (ctypes.c_double * 1)(float(weight))
            weight_arg = ctypes.cast(weight_array, dblptr)
        else:
            weight_arg = a2c(keepalive,weight)

        c_polx = a2c(keepalive,polx) if polx is not None else None
        c_poly = a2c(keepalive,poly) if poly is not None else None
        c_polz = a2c(keepalive,polz) if polz is not None else None
        c_userflags = ( a2c(keepalive,pdgcode,'uint32',uint32ptr)
                        if userflags is not None else None )
        rawfct_wm( _str2cstr(filename),
                   ctypes.c_uint64(nparticles),
                   ctypes.c_uint(flags),
                   a2c(keepalive,ekin),
                   a2c(keepalive,x), a2c(keepalive,y), a2c(keepalive,z),
                   a2c(keepalive,ux), a2c(keepalive,uy), a2c(keepalive,uz),
                   a2c(keepalive,time), weight_arg, pdgcode_arg,
                   c_polx, c_poly, c_polz, c_userflags )
    allfcts['write_mcpl'] = fct_wm

    __cache_fcts[0]=allfcts
    return __cache_fcts[0]
