
def parse_args(argv=None):
    from argparse import ArgumentParser
    import textwrap
    def wrap(t,w=59):
        return textwrap.fill( ' '.join(t.split()), width=w )

    parser = ArgumentParser(
        description=wrap("""Generate new MCPL files based upon an existing
                            kernel density estimate XML file.""",
                         79 )
    )
    parser.add_argument('kdefile',metavar='KDEFILE',
                        help=wrap("""Kernel Density Estimate .xml file generated
                        via KDSource python interface."""))
    parser.add_argument('outmcpl',metavar='OUTMCPL',
                        help=wrap("""Name of MCPL file which will be created."""))
    parser.add_argument('nparticles',metavar='NPARTICLES',type=int,
                        help=wrap("""Number of particles to sample and place in
                        OUTMCPL."""))
    parser.add_argument('-s','--seed',metavar='SEED',type=int,
                        help=wrap("""Seed value of Number of particles to sample
                        and place in OUTMCPL. If not set, a random seed is
                        generated (based on bytes from os.urandom)."""))
    parser.add_argument('-f','--force',action='store_true',
                        help=wrap("""Will overwrite existing file if it already
                        exists."""))

    return parser.parse_args(argv)

def main(argv=None):
    from .resample import resample_to_mcpl
    args = parse_args(argv)
    resample_to_mcpl( args.kdefile, args.outmcpl, args.nparticles,
                      force=args.force, rng = args.seed )

if __name__ == '__main__':
    main()
