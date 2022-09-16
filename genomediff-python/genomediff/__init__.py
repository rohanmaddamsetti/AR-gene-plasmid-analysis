import itertools
from genomediff.parser import GenomeDiffParser
from genomediff.records import Metadata


class GenomeDiff(object):

    def __init__(self):
        self.metadata = {}
        self.mutations = []
        self.evidence = []
        self.validation = []
        self._index = {}

    @classmethod
    def read(cls, fsock):
        gd = GenomeDiff()

        for record in GenomeDiffParser(document=gd, fsock=fsock):
            if isinstance(record, Metadata):
                gd.metadata[record.name] = record.value
            else:
                if len(record.type) == 3:
                    gd.mutations.append(record)
                if len(record.type) == 2:
                    gd.evidence.append(record)
                if len(record.type) == 4:
                    gd.validation.append(record)
                gd._index[record.id] = record
        return gd

    def __getitem__(self, item):
        return self._index[item]

    def write(self, fsock):
        raise NotImplementedError()

    def __len__(self):
        return len(self.mutations) + len(self.evidence) + len(self.validation)

    def __iter__(self):
        return itertools.chain(self.mutations, self.evidence, self.validation)

    def __str__(self):
        return '\n'.join(["MUTATIONS:",'\n'.join([str(x) for x in self.mutations]),
                          "EVIDENCE:",'\n'.join([str(x) for x in self.evidence]),
                          "VALIDATION:",'\n'.join(self.validation)])


    def remove(self,*args, mut_type=None):
        ''' 
        Remove mutations that satisfy the given conditions. Implementation of
        gdtools REMOVE for genomediff objects.
        
        Input: a variable number of conditions, e.g. 'gene_name==rrlA','frequency>=0.9'.
               If mut_type is specified, only that mutation type will be removed.
        Output: self.mutations is updated, with mutations satifying the conditions
                having been removed.
        '''
        updated_mutations = []
        for rec in self.mutations:
            if (mut_type is None or mut_type == rec.type) and rec.satisfies(*args):
                continue       
            else:
                updated_mutations.append(rec)

        self.mutations = updated_mutations
