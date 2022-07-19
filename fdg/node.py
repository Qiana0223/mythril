class Node():
    def __init__(self,selector:str, index:int, full_name:str):
        self.selector=selector
        self.index=index
        self.full_name=full_name


class DirectEdge():
    def __init__(self,tail: str, head:str, edgeLabel: str):
        self.tail=tail
        self.head=head
        self.edge_label=edgeLabel
