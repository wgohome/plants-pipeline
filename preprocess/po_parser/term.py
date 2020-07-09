from collections import defaultdict
import networkx as nx
import re
# local imports
from .relation import Relation

class Term():
    def __init__(self, ontology, string):
        self.ontology = ontology
        attrs = self.unpack(string)
        self.__dict__.update(**attrs)
        self.parents = [] # Direct parents only
        self.children = [] # Direct children only

    def __getattr__(self, name):
        return []

    def unpack(self, string):
        lines = string.strip().split('\n')[1:]
        attrs = defaultdict(list)
        for line in lines:
            key = re.findall(r'([\w]+?):', line)[0]
            val = re.findall(r': (.+)$', line)[0]
            attrs[key].append(val)
        return attrs

    @property
    def parent_names(self):
        return [parent.name for parent in self.parents]

    @property
    def children_names(self):
        return [children.name for children in self.children]

    def superclasses(self):
        """Returns a list of lists (parent paths of Terms).
        Only considers 'is_a' relations."""
        # Recursive implementation of finding paths
        target = self.ontology.find_by_name('Plant Anatomical Entity')
        return self.ontology.simple_paths(self, target)

    def print_superclasses(self):
        target = self.ontology.find_by_name('Plant Anatomical Entity')
        return self.ontology.print_simple_paths(self, target)

    def subclasses(self):
        full_paths = []
        paths = [[child] for child in self.children]
        while paths != []:
            for path in paths:
                # pdb.set_trace()
                paths.remove(path)
                result, leaf = self.expand_path(path)
                if leaf:
                    full_paths.append(result)
                else:
                    paths.extend(result)
        return full_paths

    def print_subclasses(self):
        for path in self.subclasses():
            print('•', ' ← '.join(term.name[0] for term in path))

    @staticmethod
    def expand_path(path):
        children = path[-1].children
        if children != []:
            return [path + [child] for child in children], False
        else:
            return path, True
