from contextlib import contextmanager

@contextmanager
def work_in_tmpdir():
    """Context manager for working in a temporary directory (automatically
    created+cleaned) and then switching back"""
    import os
    import tempfile
    the_cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(the_cwd)#Important to leave tmpdir *before* deletion, to
                             #avoid PermissionError on Windows.

def _test_samples_mcpl_gz_data():
    data = ('H4sICEXhF2EAA3NhbXBsZXMyLm1jcGwA83UO8DEwMPbhZMAOGKG0CpQGqXMMdvb0'
            'VAgODru55siBYjdjB5CEyK0jtjHd6fs+XKiyBfHZdxTZO3AwMBQ9PHqAyzfZESS2'
            '3+vb3qUdM21fBoTagfh2pRVgNaf4pjmyzNlsDxLzdZtqe+rA0r2fd96yAfHnyC0E'
            'q3ES2HCgsTJ7P0jMrWDKHqUiUTufpHVgu0K0e8BqsqJMHLeciAHbVZcevlflkdre'
            'N3KnwGr4b6eA1fycIOBwwVcbbNec8nm7md+J7i15FAB2z6uGOrCa5cU8B9t9ZcHm'
            '8GVH253r3bFnhYEi2JwHxQ1gNTZ7kh03V6wDu+ejpNnuH6e37D0asQ2spkwH4veU'
            'V5oOohW7wOGjv0LFJmD2DNscCz2wXdErZoPVAADw9C1tgQEAAA==')
    import base64
    return base64.b64decode(data)

def _test_source_xml_data():
    data = ('PD94bWwgdmVyc2lvbj0iMS4wIiA/Pgo8S0RTb3VyY2U+Cgk8SiB1bml0cz0iMS9'
            'zIj4xLjA8L0o+Cgk8UExpc3Q+CgkJPHB0Pm48L3B0PgoJCTxtY3BsbmFtZT5zYW'
            '1wbGVzLm1jcGwuZ3o8L21jcGxuYW1lPgoJCTx0cmFzbC8+CgkJPHJvdC8+CgkJP'
            'Hgyej4wPC94Mno+Cgk8L1BMaXN0PgoJPEdlb20gb3JkZXI9IjMiPgoJCTxMZXRo'
            'YXJneT4KCQkJPGRpbT4xPC9kaW0+CgkJCTxwYXJhbXMgbnBzPSIxIj4xMDwvcGF'
            'yYW1zPgoJCTwvTGV0aGFyZ3k+CgkJPFN1cmZYWT4KCQkJPGRpbT4yPC9kaW0+Cg'
            'kJCTxwYXJhbXMgbnBzPSI1Ij4taW5mIGluZiAtaW5mIGluZiAwPC9wYXJhbXM+C'
            'gkJPC9TdXJmWFk+CgkJPFBvbGFyTXU+CgkJCTxkaW0+MjwvZGltPgoJCQk8cGFy'
            'YW1zIG5wcz0iMCIvPgoJCTwvUG9sYXJNdT4KCQk8dHJhc2wvPgoJCTxyb3QvPgo'
            'JPC9HZW9tPgoJPHNjYWxpbmc+Ni45Njg2ODE0NGUtMDEgOS42NzQ0MDYzMGUrMD'
            'AgNS41OTQxNTYwMmUrMDAgOS41MjU5NjkxNWUtMDMKIDkuMzU3MDUxNDZlKzAxP'
            'C9zY2FsaW5nPgoJPEJXIHZhcmlhYmxlPSIwIj4wLjczNzIzNjQ2Njk4MTAzMTU8'
            'L0JXPgo8L0tEU291cmNlPgo=')
    import base64
    return base64.b64decode(data)

def test_resample():
    import pathlib
    import mcpl
    import numpy
    from .resample import resample_to_mcpl

    pathlib.Path('samples.mcpl.gz').write_bytes(_test_samples_mcpl_gz_data())
    pathlib.Path('source.xml').write_bytes(_test_source_xml_data())
    outfn = 'foo_tmp_testfile.mcpl.gz'
    outpath = pathlib.Path(outfn)
    if outpath.is_file():
        outpath.unlink()
    resample_to_mcpl( 'source.xml', outfn, 117, force=True )
    if not outpath.is_file():
        raise RuntimeError('Test failed: did not create resampled MCPL file')
    f = mcpl.MCPLFile('foo_tmp_testfile.mcpl.gz')
    if not f.nparticles == 117:
        if not outpath.is_file():
            raise RuntimeError('Test failed: Unexpected nparticles in'
                               ' resampled MCPL file')
    for p in f.particle_blocks:
        if not numpy.all(p.pdgcode == 2112):
            raise RuntimeError('Test failed: Unexpected pdgcode in'
                               ' resampled MCPL file')
        for a in ['z','time','weight']:
            v = 0.0 if a!='weight' else 1.0
            if not numpy.all(getattr(p,a) == v):
                raise RuntimeError(f'Test failed: Unexpected {a} in'
                                   ' resampled MCPL file')

def test_importpyapi():
    try:
        from .api import KDSource # noqa F401
        from .api import load # noqa F401
    except ImportError as e:
        raise RuntimeError('Could not import kdsource.api') from e

def test_all():
    print("==> Performing basic test of kdsource installation")
    print("  ==> Testing resampling")
    test_resample()
    print("  ==> Testing dependencies for Python API for training")
    test_importpyapi()
    print("==> All tests ended OK")

def main(argv = None):
    if argv is None:
        import sys
        argv = sys.argv
    if '--cwd' in argv[1:]:
        test_all()
    else:
        with work_in_tmpdir():
            test_all()

if __name__ == '__main__':
    main()
