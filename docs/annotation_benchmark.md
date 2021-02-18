# G. Annotation Benchmark

## Percentage Accuracy

The script to compare the annotations from LSTrAP to the ones from EVOREPRO is found [here](https://github.com/wirriamm/plants-pipeline/blob/master/preprocess/benchmark.py). The comparison is done by matching the PO term of the LSTrAP-Kingdom annotation to the EVOREPRO term, where each of the 10 possible EVOREPRO organ terms are matched to their respective PO term.

The mapping of the EVOREPRO term to PO term is given in the script as such:
```
po_mapper = {
    'Spore': ['plant spore'],
    'Female': ['plant egg cell', 'plant ovary', 'plant ovule'],
    'Leaf': ['leaf'],
    'Stem': ['stem'],
    'Male': ['pollen', 'microspore'],
    'Flower': ['flower'],
    'Seeds': ['seed'], # is a leaf node. No sublasses.
    'Root': ['root'], # does not have any meristem as children
    'Root meristem': ['root meristem'], # path does not coincide with root po term
    'Apical meristem': ['apical meristem']
}
```

This will result in a combined table, which is further edited to check for matching terms between LSTrAP-Kingdom and EVOREPRO, but is not identified by using Plant Ontology directed acyclic graph relations. The resulting edited table is provided in Supplementary Material Table S5.

## Percentage coverage

The percetage coverage of the annotation of each species is given by the number of samples annotated by LSTrAP-kingdom over the total number of samples processed, for each species.
