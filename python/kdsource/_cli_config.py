
def parse_args(argv=None):
    from argparse import ArgumentParser
    from .config import config_items
    import textwrap
    def wrap(t,w=59):
        return textwrap.fill( ' '.join(t.split()), width=w )

    parser = ArgumentParser(
        description
        ='Provide basic information about KDSource installation'
    )
    parser.add_argument('-v','--version',action='store_true',
                        help=wrap("""Show the KDSource version number and
                        exit."""))
    parser.add_argument('-i','--intversion',action='store_true',
                        help=wrap("""Show KDSource version encoded into
                        single integral number (e.g. v1.2.3 is 1002003) and
                        exit."""))
    parser.add_argument('-s','--summary',action='store_true',
                        help=wrap("""Print summary information about
                        installation and exit.  This displays all the
                        information that is otherwise available via the --show
                        flag."""))
    parser.add_argument("--show",metavar="ITEM", choices=['list']+config_items,
                        help=wrap("""Print value of the requested information
                        ITEM for the current NCrystal installation and exit. Run
                        with "--show list" to get a list of available ITEM
                        values."""))
    return parser.parse_args(argv)


def main(argv=None):
    from .config import summary, config, config_items
    args = parse_args(argv)
    if args.version:
        print(config('version'))
    elif args.intversion:
        print(config('intversion'))
    elif args.summary:
        summary()
    elif args.show:
        if args.show == 'list':
            for c in sorted(config_items):
                print(c)
        elif args.show not in config_items:
            raise SystemExit('ERROR: Unknown configuration item "%s"'%args.show)
        else:
            print(config(args.show))
    else:
        summary()


if __name__ == '__main__':
    main()
