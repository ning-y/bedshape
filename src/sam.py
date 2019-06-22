import logging
from alignment import Alignment, CigarUnavailableError, ClippedRegionEmptyError

logger = logging.getLogger('bedshape')

class File:
    """
    Represents a SAM file.
    """

    def __init__(self, iterable_lines):
        lines_skipped = 0
        total_lines = 0
        self.lines = []
        for line in iterable_lines:
            total_lines += 1
            try:
                self.lines.append(Line(line))
            except CigarUnavailableError:
                lines_skipped += 1
                logger.debug('skipped {} (no mapping information)'.format(
                        line.split()[0]))

        logger.warn('{} lines were skipped (not mapping to the region)'.format(
                lines_skipped))
        logger.info('Read {} lines ({} skipped)'.format(
                total_lines-lines_skipped, lines_skipped))

    def soft_clip(self, start, stop):
        lines_skipped = 0
        total_lines = 0
        clipped_lines = []
        for line in self.lines:
            total_lines += 1
            try:
                line.soft_clip(start, stop)
                clipped_lines.append(line)
            except ClippedRegionEmptyError:
                lines_skipped += 1
                logger.debug('skipped {} (no mappable bases after clipping)'.format(
                        line.fields[0]))
        logger.warn('{} lines did not map after the clip'.format(lines_skipped))
        logger.info('Clipped {} lines ({} skipped)'.format(
                total_lines-lines_skipped, lines_skipped))
        self.lines = clipped_lines

    def __repr__(self):
        return '\n'.join([str(line) for line in self.lines])

class Line:
    """
    Represents a line in the SAM file.
    """

    TYPE_HEADER = 0
    TYPE_ALIGNMENT = 1

    def __init__(self, line_string):
        self.type = self.TYPE_HEADER if line_string.startswith('@') \
                else self.TYPE_ALIGNMENT

        if self.type == self.TYPE_HEADER:
            self.fields = [line_string]
            return

        self.fields = line_string.split()
        pos, cigar = self.fields[3], self.fields[5]

        if cigar == '*':
            raise CigarUnavailableError

        md = next(filter(lambda field: field.startswith('MD:Z:'), self.fields))
        md = md.replace('MD:Z:', '')
        self.alignment = Alignment(pos, cigar, md)

    def soft_clip(self, start, stop):
        if self.type == self.TYPE_HEADER:
            return

        self.strip_paired_end_info()

        self.fields[2] = '{}:{}-{}'.format(self.fields[2], start, stop)
        self.alignment.soft_clip(start, stop)
        self.fields[3] = str(self.alignment.pos)
        self.fields[5] = self.alignment.cigar
        self.fields = list(map(
                lambda field: 'MD:Z:'+self.alignment.md if \
                        field.startswith('MD:Z:') else field,
                self.fields))

    def strip_paired_end_info(self):
        '''
        fields[1]: Bitwise flags according to the SAM specifications:
              1 -- template having multiple segments in sequencing
              2 -- each segment properly aligned according to the aligner
              4 -- segment unmapped
              8 -- next segment in template unmapped
             16 -- SEQ being reverse complemented
             32 -- SEQ of the next segment in the template being reversed
                   complemented
             64 -- the first segment in the template
            128 -- the last segment in the template
            ...
        fields[6]: reference sequence name of the primary alignment of the next
                   read in the template; '*' when information is unavailable.
        fields[7]: 1-based position of the primary alignment of the next read in
                   the template; '0' when information is unavailable.
        fields[8]: signed observed template length; '0' for single-segment
                   template, or when information is unavailable.
        '''
        flags = int(self.fields[1])
        flags &= 0b00111100
        self.fields[1] = str(flags)

        self.fields[6:9] = ['*', '0', '0']

    def __repr__(self):
        return '\t'.join(self.fields)


class NoMappableBaseException(Exception):
    """
    Exception raised if after a soft-clip operation, there are no more mappable
    bases in the transcript (e.g. all deletions).
    """
