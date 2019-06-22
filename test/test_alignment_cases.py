import os, sys
import pytest

sys.path.append(os.path.join(sys.path[0], '../src'))

from alignment import Alignment

bowtie2_cigars = [
    ('148M',          [('M', 148)]),
    ('23M1D95M',      [('M', 23), ('D', 1), ('M', 95)]),
    ('11M1D73M2D65M', [('M', 11), ('D', 1), ('M', 73), ('D', 2), ('M', 65)]),
    ('30M2D119M1S',   [('M', 30), ('D', 2), ('M', 119), ('S', 1)]),
    ('9S105M2D35M',   [('S', 9), ('M', 105), ('D', 2), ('M', 35)]),
    ('27M2D122M',     [('M', 27), ('D', 2), ('M', 122)]),
    ('147M1S',        [('M', 147), ('S', 1)]),
    ('150M',          [('M', 150)]),
    ('88M2D24M2D38M', [('M', 88), ('D', 2), ('M', 24), ('D', 2), ('M', 38)])
]

hisat2_cigars = [
    ('43M132N60M47S',     [('M', 43), ('N', 132), ('M', 60), ('S', 47)]),
    ('30M89N103M132N17M', [('M', 30), ('N', 89), ('M', 103), ('N', 132), ('M', 17)])
]

bowtie2_mds = [
    ('7^GC2^TG4G3C86', [
        [Alignment.MD_MATCH, 7], [Alignment.MD_DELETION, 'G', 'C'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_DELETION, 'T', 'G'],
        [Alignment.MD_MATCH, 4], [Alignment.MD_MISMATCH, 'G'],
        [Alignment.MD_MATCH, 3], [Alignment.MD_MISMATCH, 'C'],
        [Alignment.MD_MATCH, 86]]),
    ('2A2G3A2G0^A71', [
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'G'],
        [Alignment.MD_MATCH, 3], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'G'],
        [Alignment.MD_DELETION, 'A'], [Alignment.MD_MATCH, 71]]),
    ('4^TA3C96', [
        [Alignment.MD_MATCH, 4], [Alignment.MD_DELETION, 'T', 'A'],
        [Alignment.MD_MATCH, 3], [Alignment.MD_MISMATCH, 'C'],
        [Alignment.MD_MATCH, 96]]),
    ('2C0T2A1A2A1A1C0T0A0G0T5^GT4', [
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'C', 'T'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 1], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 1], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 1], [Alignment.MD_MISMATCH, 'C', 'T', 'A', 'G', 'T'],
        [Alignment.MD_MATCH, 5], [Alignment.MD_DELETION, 'G', 'T'],
        [Alignment.MD_MATCH, 4]]),
    ('2A2G3A2G0^A71', [
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'G'],
        [Alignment.MD_MATCH, 3], [Alignment.MD_MISMATCH, 'A'],
        [Alignment.MD_MATCH, 2], [Alignment.MD_MISMATCH, 'G'],
        [Alignment.MD_DELETION, 'A'], [Alignment.MD_MATCH, 71]])
]

bowtie2_fields = [
    ('1M', '1'), ('10M', '10'), ('5S10M', '10'), ('10M5S', '10'),
    ('5S10M5S', '10'), ('5S10M5S', '4C5'), ('5S10M5S', '3A0C0T0G3'),
    ('5S10M5S', 'A0C6T0G'), ('1M1D', '1^A'), ('1D1M', '^A1'), ('1M1D1M', '1^A1'),
    ('1M1D1M', '1^A0A'), ('1M1D1M', 'A0^A0A'), ('23M1D95M', '23^A85C9'),
    ('11M1D73M2D65M', '11^A73^AC65'), ('30M2D119M1S', '30^AA83G35'),
    ('9S105M2D35M', '7C0T96^AA35'), ('27M2D122M', '27^AA122'),
    ('147M1S', '21A125'), ('150M', '47C11A15C24A49'),
    ('88M2D24M2D38M', '47A15C24^AC24^AA38'), ('150M', '129A20')
]

hisat2_fields = [
    ('150M', '137T12'), ('145M5S', '43A14G86'), ('150M', '43A14G86G4'),
    ('1S149M', '122G26'), ('7S41M3D102M', '41^AGT76T25'),
    ('34M1D19M2D75M22S', '34^T0C18^TT7T67'), ('150M', '126T23'),
    ('46S90M14S', '71G18'), ('51M2D99M', '10C29G7T2^TG99'),
    ('108M1D42M', '48T59^T42'), ('31M2D119M', '31^TC119'),
    ('150M', '90T59'), ('150M', '27A122'), ('26S124M', '124'),
    ('13S129M345N1M7S', '130'), ('5S145M', '9A13T8A30C63T17'),
    ('34M5D116M', '34^TGGTC0T1G9T19G83'), ('150M', '50G68G30'),
    ('111M2D35M346N4M', '19G91^TA0T28G9'), ('132M18S', '96G35'),
    ('81M2D69M', '20T60^TG69')
]

clip_test_set = [
    [(1, '1M', '1'), (1, 1), (1, '1M', '1')],
    [(1, '2M', '2'), (1, 1), (1, '1M1S', '1')],
    [(1, '2M', '2'), (2, 2), (1, '1S1M', '1')],
    [(65505781, '105M41S', '105'),
        (65505800, 65505900), (1, '19S86M41S', '86')],
    [(65505695, '11M1D73M2D65M', '11^A73^AC65'),
        (65505800, 65505900), (1, '102S47M', '47')],
    [(65505841, '150M', '150'),
        (65505800, 65505900), (42, '60M90S', '60')]
]
