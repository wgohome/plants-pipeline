from itertools import cycle
import networkx as nx
import re
# local imports
from .term import Term

class Ontology():
    def __init__(self, string):
        self.terms = []
        self.id_hash = {} # To optimize term look-ups
        self.name_hash = {} # To optimize term look-ups
        for match in self.find_matches_in_string(string):
            term = Term(self, match)
            self.terms.append(term)
            self.id_hash[term.id[0]] = term # For faster lookup
            self.name_hash[term.name[0]] = term # For faster lookup
        self.update_neighbors() # Updates terms' instance variables
        self.graph = self.create_digraph()

    def find_matches_in_string(self, string):
        matches = re.findall(r"(\[Term\].+?\n\n)", string, flags=re.M|re.S)
        return matches

    def find(self, id):
        return self.id_hash[id]

    def find_by_name(self, name):
        return self.name_hash[name.lower()]

    def update_neighbors(self):
        # Can only be called after all terms have been initiated
        self.update_parents()
        self.update_children()

    def update_parents(self):
        get_id = lambda string: re.findall(r'PO:[0-9]+?\b', string)[0]
        for term in self.terms:
            term.parents = [self.find(get_id(string)) for string in term.is_a if ('PO' in string)]

    def update_children(self):
        for term in self.terms:
            term.children = [query for query in self.terms if term in query.parents]

    # TODO: write tsv file for cytoscape

    def create_digraph(self):
        graph = nx.DiGraph()
        # Add edges based on parents of each term, related by is_a
        for term in self.terms:
            graph.add_edges_from([*zip(cycle([term]),term.parents)], relation='is_a')
        return graph

    def simple_paths(self, source, target):
        """Returns list of all possible paths (as list of Terms) form source to target."""
        # Set target as 'Plant Anatomical Entity' to find the paths all the way to the top
        paths = nx.all_simple_paths(
            self.graph,
            source=source,
            target=target
        )
        return paths

    def print_simple_paths(self, source, target):
        for path in self.simple_paths(source, target):
            print('•', ' ➔ '.join(term.name[0] for term in path))
