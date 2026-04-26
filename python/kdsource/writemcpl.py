
def write_mcpl( filename, *, ekin, x, y, z, ux, uy, uz, time, weight, pdgcode,
                polx = None, poly = None, polz = None, userflags = None,
                double_prec = False, force=False ):
    """Create MCPL file in the path indicated, based on particle data in numpy
    arrays. The arrays must of course all have equal length, with the exception
    that weight and pdgcode parameters are allowed to have a single scalar
    value, representing a fixed value of all particules. Additionally, polx,
    poly, polz, and userflags can be None to leave such field out of the
    file. The double_prec flag can be used to enable double precision storage of
    floating point values, and the force flag allows overwriting existing files.

    Note that the extension of the produced file will always end up as .mcpl.gz,
    so it is least surprising if the provided filename already has such an
    extension.
    """
    #weights and pdgcode can be scalars => enable universal value
    import numbers
    universal_pdgcode = isinstance(pdgcode, numbers.Integral)
    universal_weight = isinstance(weight, numbers.Number)
    arrays = [ ekin, x, y, z, ux, uy, uz, time ]
    arraylens = [ len(e) for e in arrays ]
    if (int(polx is None)+int(poly is None)+int(polz is None)) not in (0,3):
        raise RuntimeError('Please provide all or none of polx, poly, and polz')
    if polx is not None:
        arraylens.append(len(polx))
    if poly is not None:
        arraylens.append(len(poly))
    if polz is not None:
        arraylens.append(len(polz))
    if userflags is not None:
        arraylens.append(len(userflags))

    if not universal_weight:
        arraylens.append(len(weight))
    if not universal_pdgcode:
        arraylens.append(len(pdgcode))
    nparticles=list(set(arraylens))
    if len(nparticles) != 1:
        raise RuntimeError('Not all provided data arrays have the same length')
    nparticles = int(nparticles[0])
    if nparticles < 0 or nparticles > 1e16:
        raise RuntimeError('Invalid (or unreasonable) number of'
                           f' particles requested: {nparticles}')
    filename = _check_new_mcpl_filename(filename,force=force)
    from ._chooks import _load_fcts
    _load_fcts()['write_mcpl']( filename, nparticles,
                                *arrays, weight, pdgcode,
                                polx, poly, polz, userflags,
                                double_prec = double_prec )

def _check_new_mcpl_filename(filename,force = False):
    import pathlib
    dst = pathlib.Path(filename)
    if dst.name.endswith('.mcpl'):
        dst = dst.parent.joinpath(dst.name+'.gz')
    elif not dst.name.endswith('.mcpl.gz'):
        raise RuntimeError('Output name does not end with .mcpl or'
                           f' .mcpl.gz: {dst.name}')
    assert dst.name.endswith('.mcpl.gz')
    dst_dotmcpl = dst.parent.joinpath(dst.name[:-3])
    if dst_dotmcpl.is_file():
        if force:
            dst_dotmcpl.unlink()
        else:
            raise RuntimeError(f'Already exist: {dst_dotmcpl}')
    if dst.is_file():
        if force:
            dst.unlink()
        else:
            raise RuntimeError(f'Already exist: {dst}')
    if not dst.parent.is_dir():
        raise RuntimeError(f'Directory not found: {dst.parent}')
    return dst.absolute()
