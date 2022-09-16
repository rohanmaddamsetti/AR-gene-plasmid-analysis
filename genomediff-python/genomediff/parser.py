from collections import OrderedDict
import re
from genomediff.records import Metadata, Record


TYPE_SPECIFIC_FIELDS = {
    'SNP': ('seq_id', 'position', 'new_seq'),
    'SUB': ('seq_id', 'position', 'size', 'new_seq'),
    'DEL': ('seq_id', 'position', 'size'),
    'INS': ('seq_id', 'position', 'new_seq'),
    'MOB': ('seq_id', 'position', 'repeat_name', 'strand', 'duplication_size'),
    'AMP': ('seq_id', 'position', 'size', 'new_copy_number'),
    'CON': ('seq_id', 'position', 'size', 'region'),
    'INV': ('seq_id', 'position', 'size'),
    'RA': ('seq_id', 'position', 'insert_position', 'ref_base', 'new_base'),
    'MC': ('seq_id', 'start', 'end', 'start_range', 'end_range'),
    'JC': ('side_1_seq_id',
           'side_1_position',
           'side_1_strand',
           'side_2_seq_id',
           'side_2_position',
           'side_2_strand',
           'overlap'),
    'CN': ('seq_id', 'start', 'end', 'copy_number'),
    'UN': ('seq_id', 'start', 'end'),
    'CURA': ('expert',),
    'FPOS': ('expert',),
    'PHYL': ('gd',),
    'TSEQ': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
    'PFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
    'RFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end', 'enzyme'),
    'PFGE': ('seq_id', 'restriction_enzyme'),
    'NOTE': ('note',),
}

class GenomeDiffParser(object):
    def __init__(self, fsock=None, document=None):
        self._document = document
        self._fsock = fsock

    @staticmethod
    def _convert_value(value):
        for type_ in (int, float):
            try:
                return type_(value)
            except ValueError:
                pass

        if value == '.' or value == '':
            value = None
        return value

    def __iter__(self):
        metadata_pattern = re.compile(r'^#=(\w+)\s+(.*)$')
        mutation_pattern = re.compile(r'^(?P<type>[A-Z]{2,4})'
                                      '\t(?P<id>\d+)'
                                      '\t((?P<parent_ids>\d+(,\s*\d+)*)|\.?)'
                                      '\t(?P<extra>.+)?$')

        for i, line in enumerate(self._fsock):
            if not line:
                continue
            elif line.startswith('#'):
                match = metadata_pattern.match(line)
                if match:
                    yield Metadata(*match.group(1, 2))
            else:
                match = mutation_pattern.match(line)

                if match:
                    type = match.group('type')
                    id = int(match.group('id'))

                    parent_ids = match.group('parent_ids')
                    if parent_ids:
                        parent_ids = [int(id) for id in parent_ids.split(',')]

                    extra = match.group('extra').split('\t')
                    extra_dct = OrderedDict()

                    for name in TYPE_SPECIFIC_FIELDS[type]:
                        value = extra.pop(0)
                        extra_dct[name] = self._convert_value(value)

                    for k, v in (e.split('=', 1) for e in extra):
                        extra_dct[k] = self._convert_value(v)

                    yield Record(type, id, self._document, parent_ids, **extra_dct)
                else:
                    raise Exception('Could not parse line #{}: {}'.format(i, line))
