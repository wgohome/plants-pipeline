"""
Exports PO terms (filtered != GO terms and ⊆ plant_anatomy)
"""

from pronto import Ontology
import inflect

# Get plant ontology
po = Ontology("https://raw.githubusercontent.com/Planteome/plant-ontology/master/po.obo")
po_terms = [term for term in po.terms() if ("GO" not in term.id) and (term.name != None) and (term.namespace == "plant_anatomy")]
# Instantiate inflect engine
infl = inflect.engine()

# helper methods
def get_term(name):
    return [p for p in po_terms if p.name == name][0]

def superclasses(term):
    return [*term.superclasses()]

def subclasses(term):
    return [*term.subclasses()]

__all__ = ['po_terms', 'infl', 'get_term']
