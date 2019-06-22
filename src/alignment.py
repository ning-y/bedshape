import logging, re

logger = logging.getLogger('bedshape')

class Alignment:
    r"""
    Represents an alignment
    """
    MD_MATCH, MD_MISMATCH, MD_DELETION = 0, 1, 2
    # Matches are represented by a sequence of digits
    RE_MD_MATCH = re.compile(r'[1-9]\d*')
    # Deletes are represented by a caret followed by more than one sequence
    # character. There may be a delimiting '0' at the end if the next token is a
    # mismatch.
    RE_MD_DELETION = re.compile(r'\^[A-Z]+0?')
    # Mismatches are represented by a series of '0'-delimited sequence
    # characters. There may be an additional delimiting '0' at he end if the
    # next token is a deletion.
    RE_MD_MISMATCH = re.compile(r'0?[A-Z](0[A-Z])*0?')

    def __init__(self, pos, cigar, md):
        self.pos = int(pos)
        self.cigar = cigar
        self.md = md
        self.bases, self.length = self._get_bases(int(pos), cigar, md)

    def soft_clip(self, start, stop):
        new_start = None

        # Soft-clips are inferred from read bases not in self.bases. So, to
        # apply a soft-clip to the front, simply slice bases outside of the
        # new start-stop region.
        for i, base in enumerate(self.bases):
            # This is a base that should be included.
            if base.ref_pos != None and base.ref_pos >= start:
                new_start = base.ref_pos
                break
        self.pos = new_start - start + 1
        self.bases = self.bases[i:]

        for i, base in enumerate(self.bases):
            # This is a base that should not be included
            if base.ref_pos != None and base.ref_pos > stop:
                break
        else:
            # The loop searches for the index that should *not* be included. If
            # the else clause is reached, that means that *all* indices should
            # be included, but i stops at the last index. To include the last
            # index, increment i by one.
            i += 1
        self.bases = self.bases[:i]

        if self.bases == []:
            raise ClippedRegionEmptyError

        # update the ref_pos of bases
        for base in self.bases:
            if base.ref_pos != None:
                base.ref_pos -= start
                base.ref_pos += 1  # 1-based

        self.cigar = self.get_cigar()
        self.md = self.get_md()

    def get_cigar(self):
        r"""
        Mapping is guaranteed to start on a 'M' or 'D'.
        """
        cigar_tokens = []

        # Find only the bases which map to the reference, and also deletions
        # which happen within a mapped region.
        for base in self.bases:
            last_token = None if cigar_tokens == [] else cigar_tokens[-1]
            last_mapped = None if last_token == None or last_token[1] == None else \
                    last_token[1] + last_token[-1] - 1
            if base.get_type() == Base.ONE_TO_ONE:
                # if needs init, is new op, or was reference skip.
                if last_token == None or \
                        last_token[0] != 'M' or \
                        base.ref_pos - last_mapped > 1:
                    cigar_tokens.append(['M', base.ref_pos, base.rd_pos, 1])
                else:
                    last_token[-1] += 1
            elif base.get_type() == Base.INSERTION:
                if last_token == None or last_token[0] != 'I':
                    cigar_tokens.append(['I', None, base.rd_pos, 1])
                else:
                    last_token[-1] += 1
            elif base.get_type() == Base.DELETION:
                if cigar_tokens == [] or \
                        last_token[0] != 'D' or \
                        base.ref_pos - last_mapped > 1:
                    cigar_tokens.append(['D', base.ref_pos, None, 1])
                else:
                    last_token[-1] += 1

        # Convert the previously found mapping regions (including deletions) to
        # CIGAR string tokens, while taking note of reference skips.
        cigar_string = []
        last_mapped_ref_pos = None
        for token in cigar_tokens:
            # Detect reference skips ('N')
            if last_mapped_ref_pos != None and token[1] != None and \
                    token[1] - last_mapped_ref_pos > 1:
                cigar_string.append('{}N'.format(
                    token[1] - last_mapped_ref_pos - 1))
            # Handle mappings and deletions
            if token[0] == 'M':
                cigar_string.append('{}M'.format(token[-1]))
                last_mapped_ref_pos = token[1] + token[-1] - 1
            elif token[0] == 'D':
                cigar_string.append('{}D'.format(token[-1]))
                last_mapped_ref_pos = token[1] + token[-1] - 1
            elif token[0] == 'I':
                cigar_string.append('{}I'.format(token[-1]))

        # 'Pad' the CIGAR string with soft clips.
        first_mapped_rd_pos = min(
                map(lambda t: t[2],
                    filter(lambda t: t[2] != None, cigar_tokens)))
        left_soft_clip_reps = first_mapped_rd_pos - 1
        last_mapped_region = max(
                filter(lambda t: t[2] != None, cigar_tokens),
                key=lambda t: t[2])
        last_mapped_rd_pos = last_mapped_region[2] + last_mapped_region[-1] - 1
        right_soft_clip_reps = self.length - last_mapped_rd_pos
        if left_soft_clip_reps != 0:
            cigar_string.insert(0, '{}S'.format(left_soft_clip_reps))
        if right_soft_clip_reps != 0:
            cigar_string.append('{}S'.format(right_soft_clip_reps))

        return ''.join(cigar_string)

    def get_md(self):
        md = []

        # Collate all MD-relevant information
        for base in self.bases:
            last_md = None if md == [] else md[-1]
            is_match = base.get_type() == Base.ONE_TO_ONE and \
                    base.rd_base == None
            is_mismatch = base.get_type() == Base.ONE_TO_ONE and \
                    base.rd_base != None
            is_deletion = base.get_type() == Base.DELETION

            if is_deletion and base.ref_base == None:
                raise RuntimeError(
                        'Marked as deletion, but no info about base identity')

            if is_match and (last_md == None or last_md[0] != '='):
                md.append(['=', 1])
            elif is_match and last_md[0] == '=':
                last_md[-1] += 1
            elif is_mismatch and (last_md == None or last_md[0] != 'X'):
                md.append(['X', base.rd_base])
            elif is_mismatch and last_md[0] == 'X':
                last_md.append(base.rd_base)
            elif is_deletion and (last_md == None or last_md[0] != 'D'):
                md.append(['D', base.ref_base])
            elif is_deletion and last_md[0] == 'D':
                last_md.append(base.ref_base)
            else:
                continue

        # Collect into tokens
        md_tokens = []

        for op in md:
            last_op = None if md_tokens == [] else \
                    '=' if isinstance(md_tokens[-1], int) else \
                    'D' if md_tokens[-1].startswith('^') else \
                    'X'

            # Each op is guaranteed different from the last.
            if op[0] == '=':
                md_tokens.append(op[1])
            elif op[0] == 'D' and (last_op == None or last_op == '='):
                md_tokens.append('^{}'.format(''.join(op[1:])))
            elif op[0] == 'D' and last_op == 'X':
                md_tokens[-1] += '0'
                md_tokens.append('^{}'.format(''.join(op[1:])))
            elif op[0] == 'X' and (last_op == None or last_op == '='):
                md_tokens.append('0'.join(op[1:]))
            elif op[0] == 'X' and last_op == 'D':
                md_tokens[-1] += '0'
                md_tokens.append('0'.join(op[1:]))
            else:
                raise RuntimeError('Error writing the MD field')

        return ''.join(map(str, md_tokens))

    def __eq__(self, other):
        if not isinstance(other, Alignment):
            return False
        if len(self.bases) != len(other.bases):
            return False
        for b1, b2 in zip(self.bases, other.bases):
            if b1 != b2:
                return False
        return self.cigar == other.cigar and \
                self.md == other.md and \
                self.pos == other.pos

    def __repr__(self):
        return 'pos = {}; cigar = {}; md = {};\n\tbases: {}'.format(
                self.pos, self.get_cigar(), self.get_md(), ', '.join([
                    str(b) for b in self.bases]))

    @classmethod
    def _get_bases(cls, pos, cigar, md):
        r"""
        Obtain a ``Base`` list representation of this alignment from its SAM
        position, CIGAR, and MD fields.

        Returns
        -------
        list of ``Base``\ s
        """
        _cigar, _md = cigar, md
        cigar = cls._parse_cigar(cigar)
        length = sum(map(lambda t: t[1],
                filter(lambda t: t[0] in 'MIS=X', cigar)))

        # Obtain information from the CIGAR field.
        rd_pos = 1
        ref_pos = pos
        bases = []
        for token in cigar:
            op, reps = token

            # Hard clips and paddings are not relevant to our purposes.
            if op in ('H', 'P'):
                continue
            elif op in ('M', '=', 'X'):
                read_span = range(rd_pos, rd_pos+reps)
                ref_span = range(ref_pos, ref_pos+reps)
                spans = zip(read_span, ref_span)
                bases.extend([
                    Base(rd_pos=s[0], ref_pos=s[1]) for s in spans])
                rd_pos += reps
                ref_pos += reps
            elif op == 'I':
                read_span = range(rd_pos, rd_pos+reps)
                bases.extend([Base(rd_pos=s) for s in read_span])
                rd_pos += reps
            elif op == 'S':
                rd_pos += reps
            # Treat deletions and skips differently, because MD treats them
            # differently---deletions are recorded in the MD field, but skips
            # are not.
            elif op == 'D':
                ref_span = range(ref_pos, ref_pos+reps)
                bases.extend([Base(ref_pos=s) for s in ref_span])
                ref_pos += reps
            elif op == 'N':
                ref_pos += reps

        # Obtain information from the MD field.
        md = cls._parse_md(md)
        for base in bases:
            # The MD field does not contain information about insertions (I) or
            # soft clips (S)
            if base.get_type() == Base.INSERTION:
                continue
            # MD contains no additional information about matches, carry on.
            elif base.get_type() == Base.ONE_TO_ONE and \
                    md[0][0] == cls.MD_MATCH:  # FIXME IndexError
                md[0][1] -= 1
            # MD contains additional information about mismatches, which is the
            # base identity of the nucleotide on the read.
            elif base.get_type() == Base.ONE_TO_ONE and \
                    md[0][0] == cls.MD_MISMATCH:
                base.rd_base = md[0][1]
                md[0].pop(1)
            # MD contains additional information about deletions, which is the
            # base identity of the nucleotide on the reference (which was
            # deleted from the read).
            elif base.get_type() == Base.DELETION and \
                    md[0][0] == cls.MD_DELETION:
                base.ref_base = md[0][1]
                md[0].pop(1)
            else:
                raise RuntimeError(
                        'CIGAR and MD fields do not match: ({}, {})'.format(
                            _cigar, _md))

            # Now, remove information already read from the MD tokens.
            if md[0][0] == cls.MD_MATCH and md[0][1] <= 0:
                md = md[1:]
            elif md[0][0] in (cls.MD_MISMATCH, cls.MD_DELETION) and \
                    len(md[0]) <= 1:
                md = md[1:]

        return bases, length

    @staticmethod
    def _parse_cigar(cigar):
        r"""
        Tokenises a CIGAR string.

        Returns
        -------
        list of two-element tuples of str (op) and int (reps)
        """
        if cigar == '*':
            raise CigarUnavailableError

        tokens = []
        reps = ''
        for char in cigar:
            if char.isnumeric():
                reps += char
            else:
                tokens.append((char, int(reps)))
                reps = ''
        assert reps == ''  # must end with non-numeric

        return tokens

    @classmethod
    def _parse_md(cls, md):
        r"""
        Tokenises a MD string.

        Returns
        -------
        list of lists of length >= 2. The first element of the list indicates
        the type of MD token, and is accessed by the static class attributes.

        If match, the list is of length 2 and the second element represents the
        number of matches. If deletion, the list is of length ``n + 1`` where n
        is the number of deletions, and the last ``n`` elements are the string
        base identities of the reference which were deleted. If mismatch, the
        list is the same as the list for deletions, except the last ``n``
        elements are the string base identities of the read which are
        mismatched.
        """
        md_tokens = []
        while md != '':
            if cls.RE_MD_MATCH.match(md):
                matching = cls.RE_MD_MATCH.match(md).group()
                md_tokens.append([cls.MD_MATCH, int(matching)])
                md = md.replace(matching, '', 1)
            elif cls.RE_MD_DELETION.match(md):
                matching = cls.RE_MD_DELETION.match(md).group()
                md_tokens.append(
                        [cls.MD_DELETION] + [
                            b for b in matching if not b in ('^', '0')])
                md = md.replace(matching, '', 1)
            elif cls.RE_MD_MISMATCH.match(md):
                matching = cls.RE_MD_MISMATCH.match(md).group()
                md_tokens.append(
                        [cls.MD_MISMATCH] + [
                            b for b in matching if not b == '0'])
                md = md.replace(matching, '', 1)
            else:
                raise RuntimeError('Invalid MD string')

        return md_tokens


class Base:
    r"""
    Represents a base of an alignment.

    Attributes
    ----------
    ONE_TO_ONE: int
    DELETION: int
    INSERTION: int
    """
    ONE_TO_ONE, DELETION, INSERTION = 0, 1, 2

    def __init__(self, rd_pos=None, ref_pos=None, rd_base=None, ref_base=None):
        assert rd_pos == None or isinstance(rd_pos, int)
        assert ref_pos == None or isinstance(ref_pos, int)
        assert rd_base == None or isinstance(rd_base, str)
        assert ref_base == None or isinstance(rd_base, str)

        self.rd_pos = rd_pos
        self.ref_pos = ref_pos
        self.rd_base = rd_base
        self.ref_base = ref_base

    def get_type(self):
        r"""
        Get the type of base. Matches and mismatches are encoded as
        ``ONE_TO_ONE``. Deletions are ``DELETION``. Both insertions and soft clips
        are represented as ``INSERTION``, since if a ``Base`` with ``rd_base`` but
        not ``ref_base`` is an insertion or a soft clip depends on its context
        in the overall alignment (if adjacent to the ends, soft clip; else
        insertion.)

        Returns
        -------
        int
        """
        if self.ref_pos != None and self.rd_pos != None:
            return self.ONE_TO_ONE
        # Insertions into the reads, which are not present in the reference. 'I'
        # in CIGAR, and unrepresented in MD.
        elif self.ref_pos == None and self.rd_pos != None:
            return self.INSERTION
        # Deletions from the reads, which are present in the reference. 'D' in
        # CIGAR, and represented by tokens starting with '^' in MD.
        elif self.ref_pos != None:
            return self.DELETION
        else:
            raise RuntimeError(
                    'Invalid internal nucleotide base representation.\n'
                    'rd_pos = {}, ref_pos = {}, rd_base = {}, ref_base = {}'
                    .format(self.rd_pos, self.ref_pos, self.rd_base,
                        self.ref_base))

    def __eq__(self, other):
        return isinstance(other, Base) and \
                self.rd_pos == other.rd_pos and \
                self.ref_pos == other.ref_pos and \
                self.rd_base == other.rd_base and \
                self.ref_base == other.ref_base

    def __str__(self):
        attrs = {}
        if self.rd_pos:
            attrs['rd_pos'] = self.rd_pos
        if self.ref_pos:
            attrs['ref_pos'] = self.ref_pos
        if self.rd_base:
            attrs['rd_base'] = self.rd_base
        if self.ref_base:
            attrs['ref_base'] = self.ref_base
        return '{} base with attrs {}'.format(
                'one-to-one' if self.get_type() == self.ONE_TO_ONE else
                'deletion' if self.get_type() == self.DELETION else
                'insertion' if self.get_type() == self.INSERTION else
                '???', attrs)

class CigarUnavailableError(Exception):
    pass

class ClippedRegionEmptyError(Exception):
    pass
