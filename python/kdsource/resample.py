
def resample_to_mcpl( kds_sourcefile,
                      destination_mcpl,
                      nout,
                      rng = None,
                      force = False ):
    """Resample nout particles from source (xml) file to a new MCPL file.

    Note that the extension of the produced file will always end up as .mcpl.gz,
    so it is least surprising if the provided destination_mcpl filename already
    has such an extension. It is an error if the output file already exists
    unless force is False.

    The rng parameter is used to control the stream of random numbers used for
    the sampling. If set to an integer value or None, it will be interpreted as
    a seed value, and a dedicated RNG stream will be created according to
    random.Random(seed_value). In case of rng=None, the seed value will be
    generated based on randomness from os.urandom. Finally, the rng parameter
    can be a function pointer to a function returning floating point values
    sampled uniformly in [0,1).
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
    #Setup RNG function:
    if rng is None:
        import os
        rng = int.from_bytes(os.urandom(8), byteorder="big")
        print( "Warning: No explicit seed provided"
               f" (choosing arbitrarily seed={rng})" )
    if isinstance( rng, int ):
        import random
        custom_stream = random.Random(rng)
        rngfct = custom_stream.random
    else:
        rngfct = rng
    #Call C implementation:
    from ._chooks import _load_fcts
    _load_fcts()['resample_to_mcpl']( rngfct,
                                      ksrc.absolute(),
                                      dst.absolute(),
                                      nout )
