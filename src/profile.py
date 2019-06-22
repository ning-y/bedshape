import logging, os, shutil, subprocess, sys, tempfile, time
import alias, constants, sam

logger = logging.getLogger('bedshape')

def run(args):
    regions = []

    if args.region:
        rname = args.region.split(':')[0]
        coords = args.region.split(':')[1]
        coords = coords.replace(',', '')
        regions.append(':'.join([rname, coords]))
    elif args.bed:
        regions = parse_bedfile(args.bed)
    else:
        # should have been caught in validate(args, print_help)
        raise RuntimeError('Invalid options for profile')
    if args.region and args.bed:
        logger.warn('Both region and BED file specified. Defaulting to region')

    if args.alias:
        lookup = alias.load_json()
        reference = lookup[args.alias]['reference']
        modified = lookup[args.alias]['modified']
        unmodified = lookup[args.alias].get('unmodified', None)
        denatured = lookup[args.alias].get('denatured', None)
    elif args.reference and args.modified:
        reference = args.reference
        modified = args.modified
        unmodified = args.unmodified
        denatured = args.denatured
    else:
        # should have been caight in validate(args, print_help)
        raise RuntimeError('Invalid options for profile')
    if args.alias and \
            (args.reference or args.modified or
                args.unmodified or args.denatured):
        logger.warn(
                'Both an alias and reference/sample filepaths specified, '
                'defaulting to alias.')

    outdir = get_abs_join(
            './' if len(regions) < 2 else time.strftime('bedshape-%Y%m%d-%H%M%s'),
            args.outdir)
    for region in regions:
        logger.info(
                'Running profile for {} with\n'
                '\tmodified: {}\n'
                '\tunmodified: {}\n'
                '\tdenatured: {}'.format(region, modified, unmodified, denatured))
        profile_once(
                reference, modified, unmodified, denatured, region,
                keep=args.keep, outdir=outdir,
                min_depth=args.min_depth, max_bg=args.max_bg,
                skip_plot=args.skip_plot, skip_shape=args.skip_shape,
                min_mapq=args.min_mapq)

def validate(args, print_help):
    to_exit = False

    if not args.region and not args.bed:
        logger.error('You must specify either a region or BED file')
        to_exit = True

    if args.alias and not to_exit:
        return  # has alias, region/bed is ok, just continue
    if not args.alias and not args.reference:
        logger.error(
                'You must specify a reference through --reference, or --alias')
        to_exit=True
    if not args.alias and not args.modified:
        logger.error(
                'You must specify a modified sample through --modified, or '
                '--alias')
        to_exit=True

    if to_exit:
        print_help()
        sys.exit(1)

def parse_bedfile(bedfile):
    regions = []
    with open(bedfile) as inputfile:
        for line in inputfile:
            if line.startswith('browser') or line.startswith('track'):
                continue
            line = line.split()
            regions.append('{}:{}-{}'.format(*line[:3]))
    return regions

def profile_once(
        reference, modified, unmodified, denatured, region, *, keep, outdir,
        min_depth, max_bg, skip_plot, skip_shape, min_mapq):

    start, stop = region.split(':')[1].split('-')
    start, stop = int(start), int(stop)
    tmpdir = tempfile.mkdtemp()

    ref_name = extract_from_reference(
            reference, region, tmpdir=tmpdir)
    modified_name = extract_from_alignment(
            modified, region, out_name='modified.sam', tmpdir=tmpdir)
    unmodified_name = extract_from_alignment(
            unmodified, region, out_name='untreated.sam', tmpdir=tmpdir)
    denatured_name = extract_from_alignment(
            denatured, region, out_name='denatured.sam', tmpdir=tmpdir)

    modified_counts = make_counts(
            modified_name, region, min_mapq=min_mapq, tmpdir=tmpdir)
    unmodified_counts = make_counts(
            unmodified_name, region, min_mapq=min_mapq, tmpdir=tmpdir)
    denatured_counts = make_counts(
            denatured_name, region, min_mapq=min_mapq, tmpdir=tmpdir)

    # if denatured not supplied, denatured_counts == None
    samples = list(filter(
            lambda fn: fn != None,
                [modified_counts, unmodified_counts, denatured_counts]))

    profile_filename = '{}.profile'.format(region)
    profile_filename = make_profile(samples, ref_name, profile_filename,
            tmpdir=tmpdir, outdir=outdir, min_depth=min_depth, max_bg=max_bg)

    if not skip_plot:
        figure_filename = '{}.pdf'.format(region)
        render_figure(profile_filename, figure_filename,
                outdir=outdir, min_depth=min_depth, max_bg=max_bg)

    if not skip_shape:
        make_shape(profile_filename, region, outdir=outdir)

    if not keep:
        shutil.rmtree(tmpdir)

def extract_from_reference(reference, region, tmpdir='./'):
    out_name = '{}.fa'.format(region.replace(':', '-'))
    out_name = get_abs_join(tmpdir, out_name)
    cmd = ['samtools', 'faidx', reference, region]
    with open(out_name, 'w') as outfile:
        logger.info(' '.join(cmd))
        subprocess.run(cmd, stdout=outfile)
    return out_name

def extract_from_alignment(alignment, region, tmpdir='./', out_name=None):
    if alignment == None:
        return None

    out_name = out_name if out_name != None else \
            '{}-{}.sam'.format(get_basename(alignment), region)
    out_name = get_abs_join(tmpdir, out_name)
    cmd = ['samtools', 'view', alignment, region]
    with open(out_name, 'w') as outfile:
        logger.info(' '.join(cmd))
        subprocess.run(cmd, stdout=outfile)

    return out_name

def make_counts(alignment, region, *, min_mapq, tmpdir='./',
        clipped_name=None, mut_name=None, out_name=None):
    if alignment == None:
        return None

    start, stop = region.split(':')[1].split('-')
    start, stop = int(start), int(stop)

    with open(alignment) as infile:
        sam_file = sam.File(infile)
    sam_file.soft_clip(start, stop)

    clipped_name = clipped_name if clipped_name != None else \
            '{}.clipped.sam'.format(get_basename(alignment))
    clipped_name = get_abs_join(tmpdir, clipped_name)
    with open(clipped_name, 'w') as outfile:
        outfile.write(str(sam_file))

    mut_name = mut_name if mut_name != None else \
            '{}.mut'.format(get_basename(alignment))
    mut_name = get_abs_join(tmpdir, mut_name)
    cmd = [constants.MUT_PARSER_BIN, '-i', clipped_name, '-o', mut_name,
            '--min_mapq', str(min_mapq), '--min_qual', str(min_mapq)]
    logger.info(' '.join(cmd))
    subprocess.run(cmd)

    out_name = out_name if out_name != None else \
            '{}.counts'.format(get_basename(alignment))
    out_name = get_abs_join(tmpdir, out_name)
    cmd = [constants.MUT_COUNTER_BIN, '-i', mut_name,
            '-c', out_name, '-w']
    logger.info(' '.join(cmd))
    subprocess.run(cmd)

    return out_name

def make_profile(samples, ref_name, out_name, *, tmpdir='./',
        outdir, min_depth, max_bg):
    intermediate_name = '{}.profile'.format(get_basename(out_name))
    intermediate_name = get_abs_join(tmpdir, intermediate_name)
    cmd = ['python3', constants.PROFILER_BIN, '--fa', ref_name, '--counts'] \
            + samples + ['--out', intermediate_name, '--mindepth',
            str(min_depth), '--maxbg', str(max_bg)]
    logger.info(' '.join(cmd))
    subprocess.run(cmd)

    out_name = get_abs_join(outdir, out_name)
    cmd = ['python3', constants.NORMALISER_BIN, '--tonorm', intermediate_name,
            '--normout', out_name]
    logger.info(' '.join(cmd))
    subprocess.run(cmd)

    return out_name

def render_figure(profile, out_name, *, outdir, min_depth, max_bg):
    out_name = get_abs_join(outdir, out_name)
    cmd = ['python3', constants.RENDERER_BIN, '--infile', profile,
            '--plot', out_name, '--mindepth', str(min_depth), '--maxbg',
            str(max_bg)]
    logger.info(' '.join(cmd))
    subprocess.run(cmd)

    return out_name

def make_shape(profile, basename, *, outdir):
    shape_filename = get_abs_join(outdir, basename+'.shape')
    map_filename = get_abs_join(outdir, basename+'.map')
    varna_filename = get_abs_join(outdir, basename+'.varna.txt')
    ribosketch_filename = get_abs_join(outdir, basename+'.ribosketch.txt')
    cmd = ['python3', constants.TAB2SHAPE_BIN, '--infile', profile,
            '--map', map_filename, '--shape', shape_filename, '--varna',
            varna_filename, '--ribosketch', ribosketch_filename]
    logger.info(' '.join(cmd))
    subprocess.run(cmd )

def get_abs_join(_dir, _fn):
    return os.path.abspath(os.path.join(_dir, _fn))

def get_basename(filename):
    return os.path.splitext(os.path.basename(filename))[0]
