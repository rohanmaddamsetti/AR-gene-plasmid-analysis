import re

class Metadata(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __repr__(self):
        return "Metadata({}, {})".format(repr(self.name), repr(self.value))

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class Record(object):
    def __init__(self, type, id, document=None, parent_ids=None, **attributes):
        self.document = document
        self.type = type
        self.id = id
        self.parent_ids = parent_ids
        self.attributes = attributes

    @property
    def parents(self):
        if not self.parent_ids is None:
            return [self.document[pid] for pid in self.parent_ids]
        else:
            return []

    def __getattr__(self, item):
        try:
            return self.attributes[item]
        except KeyError:
            raise AttributeError


    def __repr__(self):
        return "Record('{}', {}, {}, {})".format(self.type,
                                             self.id,
                                             self.parent_ids,
                                             ', '.join('{}={}'.format(k, repr(v)) for k, v in self.attributes.items()))

    def __str__(self):
        return self.__repr__()

    
    def __eq__(self, other):
        ''' this definition allows identical mutations in different genome diffs
            to be equal.'''
        return self.type == other.type and self.attributes == other.attributes

    def __ne__(self, other):
        return not self.__eq__(other)

    def satisfies(self, *args):
        '''
        Input: a variable number of conditions, e.g. 'gene_name==rrlA','frequency>=0.9'.
        Output: return true if all conditions are true (i.e. correspond to key-values in attributes.

        Find a condition that evaluates to false, otherwise return True.
        '''

        ## helper function to check if values are numbers
        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
        
        for c in args:
            assert type(c) == str, "error: supplied condition is not a string."
            condition_pattern = re.compile(r'^(?P<key>[_a-z]+)'
                                            '(?P<comp>==|!=|<|<=|>|>=)'
                                            '(?P<val>[-_a-zA-Z0-9\.]+)')
            condition_match = condition_pattern.match(c)
            assert condition_match, "the supplied condition\n"+c+"\n could not be parsed."
            cond_key = condition_match.group('key')
            cond_comp = condition_match.group('comp')
            cond_val = condition_match.group('val')

            try: ## in case the given condition is not in the attributes.
                attribute_val = self.attributes[cond_key]
            except:
                continue
            
            ## add quote marks around strings before eval. can leave numbers alone.
            if not is_number(cond_val):
                cond_val = "\'"+cond_val+"\'"

            if not is_number(attribute_val):
                attribute_val = "\'"+attribute_val+"\'"
            else: ## attribute_val is a number in this record-- convert to str for eval.
                attribute_val = str(attribute_val)
            expr = attribute_val+cond_comp+cond_val
            if not eval(expr):
                return False
        return True
