import getpass, json, logging, os, pprint

ALIAS_FILENAME = os.path.abspath(os.path.join(
        os.path.realpath(__file__), '../../alias.json'))
USER = getpass.getuser()

logger = logging.getLogger('bedshape')

def list():
    pprint.pprint(load_json())

def add(args):
    if not args.modified.startswith('/') or not args.reference.startswith('/') \
            or (args.unmodified != None and not args.unmodified.startswith('/')) \
            or (args.denatured != None and not args.denatured.startswith('/')):
        logger.warn(
                'One of the alias paths are relative, not absolute.\n\tThis can '
                'cause an error if you run bedshape from a different directory.')

    new_json = load_json()
    new_json[args.name] = {}
    new_json[args.name]['user'] = USER
    new_json[args.name]['modified'] = args.modified
    new_json[args.name]['reference'] = args.reference
    if args.unmodified:
        new_json[args.name]['unmodified'] = args.unmodified
    if args.denatured:
        new_json[args.name]['denatured'] = args.denatured

    save_json(new_json)

def rm(args):
    new_json = load_json()
    del new_json[args.name]

    save_json(new_json)

def load_json():
    if not os.path.isfile(ALIAS_FILENAME):
        save_json({})
        return {}
    with open(ALIAS_FILENAME) as alias_file:
        return json.load(alias_file)

def save_json(new_json):
    with open(ALIAS_FILENAME, 'w') as alias_file:
        json.dump(new_json, alias_file)
    os.chmod(ALIAS_FILENAME, 0o666)
