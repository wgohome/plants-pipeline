class Relation():
    all = []

    def __init__(self, relation_type, parent, child):
        self.relation_type = relation_type
        self.parent = parent
        self.child = child
        self.all.append(self)
