import logging, sys
import alias, cli, constants, profile

logging.basicConfig(format='[%(asctime)s/%(levelname)s] %(message)s')
logger = logging.getLogger('bedshape')
logger.setLevel(logging.INFO)

if __name__ == '__main__':
    root_parser, alias_parser, profile_parser = cli.get_parser()
    args = root_parser.parse_args()

    if args.which_subcommand == None:
        root_parser.print_help()
    elif args.which_subcommand == 'alias':
        alias_parser.print_help()
    elif args.which_subcommand == 'list':
        alias.list()
    elif args.which_subcommand == 'add':
        alias.add(args)
    elif args.which_subcommand == 'rm':
        alias.rm(args)
    elif args.which_subcommand == 'profile':
        profile.validate(args, profile_parser.print_help)
        profile.run(args)
    else:
        raise RuntimeError('Unrecognised subcommand.')
