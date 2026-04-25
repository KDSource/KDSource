
def resample_to_mcpl( kds_sourcefile, destination_mcpl, nout, force = False ):
    """Resample n=nout particles from source (xml) file to a new MCPL file.

    Note that the extension of the produced file will always end up as .mcpl.gz,
    so it is least surprising if the provided destination_mcpl filename already
    has such an extension.

    """
    import pathlib
    from .writemcpl import _check_new_mcpl_filename
    ksrc = pathlib.Path(kds_sourcefile)
    nout = int(nout)
    if not ksrc.is_file():
        raise RuntimeError(f'Missing file: {kds_sourcefile}')
    dst = _check_new_mcpl_filename(destination_mcpl,force=force)
    if nout <= 0 or nout > 1e16:
        raise RuntimeError('Invalid (or unreasonable) number of'
                           f' particles requested: {nout}')
    from ._chooks import _load_fcts
    _load_fcts()['resample_to_mcpl']( ksrc.absolute(), dst.absolute(), nout )
