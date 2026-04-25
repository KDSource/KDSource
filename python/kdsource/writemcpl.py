
def write_mcpl( filename, ekin, x, y, z, ux, uy, uz, time, weight, pdgcode,
                double_prec = False, force=False ):
    """FIXME"""
    #weights and pdgcode can be scalars => enable universal value
    import numbers
    universal_pdgcode = isinstance(pdgcode, numbers.Integral)
    universal_weight = isinstance(weight, numbers.Number)
    arrays = [ ekin, x, y, z, ux, uy, uz, time ]
    arraylens = [ len(e) for e in arrays ]
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
