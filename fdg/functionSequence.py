from fdg.node import Node
import fdg.utils as utils
import fdg.FDG_global
SEQUENCE_NUM_LIMIT=fdg.FDG_global.seq_num_limit # default 5
sub_SEQUENCE_NUM_LIMIT=3 # default 3,

class FunctionSequence():
    def __init__(self, ftn_node:Node=None,ftn_parents_list:list=None, sv_info:dict=None,state_change_sequences=None,no_state_change_sequences=None ):
        self.ftn_node=ftn_node
        self.ftn_parents_list=ftn_parents_list

        self.sv_info = sv_info

        self.state_change_sequences=state_change_sequences
        # the sequences ending with self_node that do not change states
        self.no_state_change_sequences=no_state_change_sequences

        self.generated_sequences={}

    def generate_sequences(self):
        # generate sequences
        self._get_sequences()
        sequence_selectors=[]
        sequences=[]
        length=0
        for key,value in self.generated_sequences.items():
            if isinstance(value,list):
                sorted_v=sorted(value,key=len)
                length+=len(sorted_v)
                if length>=SEQUENCE_NUM_LIMIT:
                    sequences+=sorted[0:(SEQUENCE_NUM_LIMIT-len(sequences))]
                    break
                else:
                    sequences+=sorted_v
        for seq in sequences:
            sequence_selectors.append([fdg.FDG_global.ftn_to_selector[ftn] for ftn in seq])
        return sequence_selectors





    def _get_sequences(self):
        label_parents = {}
        for node in self.ftn_parents_list:
            if node.edge_label not in label_parents.keys():
                label_parents[node.edge_label] = [node.full_name]
            else:
                label_parents[node.edge_label] += [node.full_name]

        num_unique_edge_labels=len(label_parents.keys())
        if (num_unique_edge_labels)>1:
            prt_shortests=self._get_parent_shortest_sequences()
        for num in range(num_unique_edge_labels):
            if num==0:
                self._get_sequences_signle_parent()
            else:
                sequences=self._get_sequences_multiple_parents(num+1,label_parents,prt_shortests)
                self.generated_sequences[str(num+1)+"_parents"]=sequences

    def _get_sequences_signle_parent(self)->list:
        def get_sequence_a_parent(parent_node:Node):
            sequences = []
            # create sequences only involving a single parent
            for state_change_sequence in self.state_change_sequences[parent_node.full_name]:
                seq = state_change_sequence + [self.ftn_node.full_name]
                if seq not in self.no_state_change_sequences:
                    if seq not in sequences:
                        sequences.append(seq)
            return sequences

        for node in self.ftn_parents_list:
            sequences = get_sequence_a_parent(node)
            if "1_parent" not in self.generated_sequences:
                self.generated_sequences['1_parent'] = sequences
            else:
                for seq in sequences:
                    if seq not in self.generated_sequences['1_parent']:
                        self.generated_sequences['1_parent'].append(seq)



    def _get_sequences_multiple_parents(self, num_parents:int,labels_parents:dict,prt_shortests:dict)->list:
        sequences=[]
        labels_combs=utils.get_combination_for_a_list(labels_parents.keys(),num_parents)
        if len(labels_combs)>sub_SEQUENCE_NUM_LIMIT:
            labels_combs_selected=self._get_edge_label_combinationns(labels_combs)
        else:
            labels_combs_selected=labels_combs
        for edge_label_comb in labels_combs_selected:
            # get the parents
            parents=[labels_parents[label] for label in edge_label_comb]
            parent_combinations=utils.get_combination(parents,num_parents)
            for prt_comb in parent_combinations:
                # get parent sequence combinations
                parent_sequences=[prt_shortests[parent] for parent in prt_comb]
                parent_sequence_combs=utils.get_combination(parent_sequences,num_parents)
                for prt_seq_comb in parent_sequence_combs:
                    prt_idx_seq_comb=[]
                    # use numbers to replace each function in the sequences
                    for seq in prt_seq_comb:
                        prt_idx_seq_comb+=[fdg.FDG_global.ftn_to_idx[ftn] for ftn in seq]
                    seq=self._get_a_topological_sequence(prt_idx_seq_comb)
                    if seq not in sequences:
                        sequences.append(seq)

        # replace the numbers back to function names
        final_sequences=[]
        for seq in sequences:
            ftn_seq=[fdg.FDG_global.idx_to_ftn[idx] for idx in seq if idx!=0] # do not consider constructor indexed 0
            ftn_seq+=self.ftn_node.full_name
            final_sequences.append(ftn_seq)
        return final_sequences

    def _get_edge_label_combinationns(self,edge_label_combinations)->list:
        weights= {}
        # get weight for each comb
        for i in range(len(edge_label_combinations)) :
            weight=0
            for label in edge_label_combinations[i]:
                weight+=self.sv_info[label][1]
            if weight not in weights.keys():
                weights[weight]=[i]
            else:
                weights[weight] += [i]
        # sort the collected wegiths
        sort_weights=sorted(weights.keys())

        # find the first several combs that have lower weight values
        combs_index=[]
        length=0
        for i in range(len(sort_weights)):
            values=weights[sort_weights[i]]
            length+=len(values)
            if length>=sub_SEQUENCE_NUM_LIMIT:
                combs_index+=values[0:(sub_SEQUENCE_NUM_LIMIT-len(combs_index))]
                break
            else:
                combs_index+=values

        return [edge_label_combinations[index] for index in range(len(edge_label_combinations)) if index in combs_index]


    def _get_parent_shortest_sequences(self)->dict:
        prt_shortests={}
        for node in self.ftn_parents_list:
            shortest_seq=[]
            sequences=self.state_change_sequences[node.full_name]
            len_seq=[len(seq) for seq in sequences]
            min_len=min(len_seq)
            for seq in sequences:
                if len(seq)==min_len:
                    shortest_seq.append(seq)
            prt_shortests[node.full_name]=shortest_seq
        return prt_shortests

    def _get_a_topological_sequence(self, sequences:list)->list:
        """
        get a topological sequence from multiple sequences
        start with the constructor
        end with the target function (ftn_idx)

        based on function indices
        :param sequences:
        :return: a sequence, each element is an index
        """
        def get_graph(sequences: list)->dict:
            """
            build a graph starting with the constructor and ending with self.ftn_idx
            :param sequences:
            :return:
            """
            graph = {}
            graph[0] = []
            graph[self.ftn_node.index]=[]
            for seq in sequences:
                if len(seq)==0:continue
                if len(seq) == 1:
                    if seq[0] not in graph[0]:
                        graph[0].append(seq[0])
                    if seq[0] not in graph.keys():
                        graph[seq[0]] = [self.ftn_idx]
                    else:
                        if self.ftn_node.index not in graph[seq[0]]:
                            graph[seq[0]] += [self.ftn_node.index]
                else:
                    ftn_start = seq[0]
                    # add the edge between the constructor to the first function in the sequence
                    if ftn_start not in graph[0]:
                        graph[0].append(ftn_start)

                    # add the edge based on the sequence, each consecutive two functions has an edge
                    for ftn in seq[1:]:
                        if ftn_start not in graph.keys():
                            graph[ftn_start] = [ftn]
                        else:
                            if ftn not in graph[ftn_start]:
                                graph[ftn_start] += [ftn]
                        ftn_start = ftn

                    # add the edge between the last function in the sequence and the target function
                    assert (ftn_start == seq[-1])
                    if ftn_start not in graph.keys():
                        graph[ftn_start] = [self.ftn_node.index]
                    else:
                        if self.ftn_node.index not in graph[ftn_start]:
                            graph[ftn_start] += [self.ftn_node.index]
            return graph

        # A recursive function used by topologicalSort
        def topologicalSortUtil(graph:dict,v,visited, stack):
            # Mark the current node as visited.
            visited[v] = True
            # Recur for all the vertices adjacent to this vertex
            for i in graph[v]:
                if visited[i] == False:
                    topologicalSortUtil(graph,i, visited, stack)
            # Push current vertex to stack which stores result
            stack.insert(0, v)

        # compute the number of nodes
        all_nodes=[]
        for seq in sequences:
            all_nodes+=seq
        num_nodes=2+len(set(all_nodes))

        # build the graph
        graph=get_graph(sequences)

        # get the path
        # Mark all the vertices as not visited
        visited = [False] * num_nodes
        stack = []
        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(num_nodes):
            if visited[i] == False:
                topologicalSortUtil(graph,i, visited, stack)


        return stack


if __name__=='__main__':
    # sequences=[[1,2,3,4],[1,2,4,3],[5,6,3],[1,3]]
    # ftnSeq=FunctionSequence(7)
    # path=ftnSeq._get_a_topological_sequence(sequences)
    # print(f'path={path}')
    print([1, 2, 4, 3][0:2])
    a=[[1,2,3],[1],[1,2]]
    a.sort(key=len)
    print(a)







