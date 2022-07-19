import time

import fdg
from fdg import utils
from fdg.FDG import FDG
from fdg.sequenceAndState import SequenceAndState
import numpy as np

class SequenceGeneration():
    def __init__(self,fdg:FDG=None, sequenceAndState:SequenceAndState=None, phase1_depth=2):
        self.FDG=fdg
        self.sequenceAndState=sequenceAndState
        self.phase1_depth=phase1_depth

    def generate_sequences(self,ftn_idx)->list:

        print(f'===================================')
        seconds_start3 = time.time()
        sequences_heuristic = self.generate_sequences_by_heuristic_method(ftn_idx)
        seconds_end3 = time.time()
        sequences_heuristic.sort(key=len)
        print(f'sequences_heuristic({ftn_idx})={sequences_heuristic}')
        print(f'time used to generate sequences_heuristic:{seconds_end3 - seconds_start3}')

        print(f'===================================')
        seconds_start4 = time.time()
        sequences_paper = self.generate_sequences_paper(ftn_idx)
        seconds_end4 = time.time()
        sequences_paper.sort(key=len)
        print(f'sequences_paper({ftn_idx})={sequences_paper}')
        print(f'time used to generate sequences_paper:{seconds_end4 - seconds_start4}')



        print(f'===================================')
        seconds_start1 = time.time()
        sequences_random_parentGroups = self.generate_sequences_by_randomly_parentGroups_parents(ftn_idx)
        seconds_end1 = time.time()
        sequences_random_parentGroups.sort(key=len)
        print(f'sequences_random_parentGroups({ftn_idx})={sequences_random_parentGroups}')
        print(f'time used to generate sequences_random_parentGroups:{seconds_end1 - seconds_start1}')

        print(f'===================================')
        seconds_start2 = time.time()
        sequences_random_parents = self.generate_sequences_by_randomly_parents(ftn_idx)
        seconds_end2 = time.time()
        sequences_random_parents.sort(key=len)
        print(f'sequences_random_parents ({ftn_idx})={sequences_random_parents }')
        print(f'time used to generate sequences_random_parents :{seconds_end2 - seconds_start2}')
        return sequences_paper


    def generate_sequences_by_heuristic_method(self,ftn_idx):
        sv_parents = self.FDG.get_parents(ftn_idx)
        if len(sv_parents)==0: return []
        parent_list = self.FDG.get_parent_list(ftn_idx)
        if len(parent_list)==0: return []
        print(f'---------------------------------')
        print(f'deep function: {self.FDG.contractInfo.get_name_from_index(ftn_idx)},{ftn_idx}')
        for sv,parent in sv_parents.items():
            print(f'sv {sv}: parent(s){parent}')


        num_SVs=len(sv_parents) # number of SVs that ftn_idx reads in conditions
        if num_SVs==1: # when only one state variable is read
            if isinstance(sv_parents,dict):
                sequences=self.parent_sequences_write_one_SV(ftn_idx,parent_list,fdg.FDG_global.seq_num_limit,self.phase1_depth)
                return sequences
            else: assert False,"data type error"
        elif num_SVs==2:
            if isinstance(sv_parents, dict):
                sv_list=list(sv_parents.keys())
                # get at least 2 sequences writing one SV
                sequences_1_SV = self.parent_sequences_write_one_SV(ftn_idx, parent_list, fdg.FDG_global.num_sequences_write_1_SV, self.phase1_depth)

                # get the sequences writing two SVs
                n=fdg.FDG_global.seq_num_limit-len(sequences_1_SV)
                sequences_2_SVs=self.merged_parent_sequences_write_two_SVs(ftn_idx,sv_list,sv_parents,n)
                for seq in sequences_1_SV:
                    if seq not in sequences_2_SVs:
                        sequences_2_SVs.append(seq)
                return sequences_2_SVs


            else:
                assert False, "data type error"
        else:
            if isinstance(sv_parents, dict):
                sv_list = list(sv_parents.keys())
                # get at least 2 sequences writing one SV
                sequences_1_SV = self.parent_sequences_write_one_SV(ftn_idx, parent_list,
                                                                    fdg.FDG_global.num_sequences_write_1_SV,
                                                                    self.phase1_depth)
                # get the sequences writing two SVs
                n = fdg.FDG_global.seq_num_limit - len(sequences_1_SV)-1
                sv_list=list(sv_parents.keys())
                sv_pairs=self._get_SV_pairs(sv_list)
                sequences_2_SVs=[]
                count=0
                for idx in range(len(sv_pairs)):
                    if count==n:break
                    sv_pair=sv_pairs[idx]
                    sequences = self.merged_parent_sequences_write_two_SVs(ftn_idx, sv_pair, sv_parents, 1)
                    if len(sequences)>0:
                        count+=1
                        sequences_2_SVs.append(sequences[0])

                # get one sequence writing all SVs
                sequences_3_or_more_SVs=self.merged_parent_sequences_write_more_SVs(ftn_idx,sv_list,sv_parents,1)
                for seq in sequences_1_SV:
                    if seq not in sequences_2_SVs:
                        sequences_2_SVs.append(seq)
                return sequences_2_SVs+ sequences_3_or_more_SVs
            else:
                assert False, "data type error"


    def parent_sequences_write_one_SV(self, ftn_idx:int,parent_list:list,n:int,min_length:int)->list:
        """
        get at most n parent sequences that write one SV and have length >= min_length
        :param parent_list: the parents to be considered
        :param n: specify the number of sequences returned
        :param min_length: the length of each sequence >= min_length
        :return:
        """

        parent_has_state_changing_sequences=[ prt for prt in parent_list if self.sequenceAndState.has_state_changing_sequences_length(prt,min_length)]
        if n==-1:
            # get all state-changing sequences with length >= min_length
            sequences = self.sequenceAndState.get_all_state_changing_sequeces_length(
                parent_has_state_changing_sequences, self.phase1_depth)
            return [seq + [ftn_idx] for seq in sequences]
        if len(parent_has_state_changing_sequences)>n:
            # randomly select n
            sequences=self.sequenceAndState.get_shortest_state_changing_sequeces_length(parent_has_state_changing_sequences,self.phase1_depth)
            sequences_selected=utils.random_select(sequences, n)
            return [seq+[ftn_idx] for seq in sequences_selected]
        else:
            # get all state-changing sequences with length >= min_length
            sequences=self.sequenceAndState.get_all_state_changing_sequeces_length(parent_has_state_changing_sequences,self.phase1_depth)

            sequences_selected = utils.random_select(sequences, n)
            return [seq + [ftn_idx] for seq in sequences_selected]


    def merged_parent_sequences_write_two_SVs(self, ftn_idx:int,sv_list:list,sv_parents:dict,n:int)->list:
        """
        return n parent sequences that write two state varibles
        :param sv_list:
        :param sv_parents:
        :param n:
        :return:
        """
        assert len(sv_list)==2
        sv_types=[self.FDG.contractInfo.get_sv_type(sv) for sv in sv_list]
        sv_type_primitive=[True for type in sv_types if type.__eq__(fdg.FDG_global.primitive_index)]
        if n==1 or len(sv_type_primitive)==0:
            sv1_seq=self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(sv_parents[sv_list[0]],1)
            sv2_seq = self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(sv_parents[sv_list[1]], 1)
        else:
            sv1_seq=[]
            sv2_seq=[]
            if sv_types[0]==fdg.FDG_global.primitive_index:
                sv1_seq=self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(sv_parents[sv_list[0]],2)
            else:
                sv1_seq=self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(sv_parents[sv_list[0]],1)

            if sv_types[1] == fdg.FDG_global.primitive_index:
                sv2_seq = self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(
                    sv_parents[sv_list[1]], 2)
            else:
                sv2_seq = self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(
                    sv_parents[sv_list[1]], 1)

        if len(sv1_seq) == 0 or len(sv2_seq) == 0:
            return []
        else:
            sv1_seq_indices=range(len(sv1_seq))
            sv2_seq_indices = range(len(sv2_seq))
            comb_indices=utils.get_combination([sv1_seq_indices,sv2_seq_indices],2)
            final_sequences=[]
            for comb in comb_indices:
                comb=list(comb) # there are two indices in the comb
                sequences=[sv1_seq[comb[0]]]+[sv2_seq[comb[1]]]
                merged_sequence=self._get_a_topological_sequence(ftn_idx,sequences)
                if len(merged_sequence)>=self.phase1_depth+1: #
                    final_sequences.append(merged_sequence)
            return utils.random_select(final_sequences,n)


    def merged_parent_sequences_write_more_SVs(self,ftn_idx:int,sv_list:list,sv_parents:dict,n:int)->list:
        """
        return n parent sequences writing 3 or more state variables.
        :param sv_list:
        :param sv_parents:
        :param n:
        :return:
        """
        assert len(sv_list) >=3
        assert n==1 # we only consider 1 sequence that write 3 or more state variables
        sequences=[]
        for sv in sv_list:
            parents=sv_parents[sv]
            seq=self.sequenceAndState.get_n_state_changing_sequences_from_multiple_parents(parents,1)
            if len(seq)>0:
                sequences.append(seq[0])
        if len(sequences)<3:# must have 3 or more sequences to be merged
            return []
        else:
            return [self._get_a_topological_sequence(ftn_idx,sequences)]


    def _get_a_topological_sequence(self,ftn_idx:int, sequences:list)->list:
        """
        get a topological sequence from multiple sequences;
        start with the constructor;
        end with the target function (ftn_idx);

        based on function indices
        :param sequences:
        :return: a sequence, each element is an index
        """
        def get_graph(sequences: list,ftn_idx:int)->dict:
            """
            build a graph starting with the constructor and ending with self.ftn_idx
            :param sequences:
            :return:
            """
            graph = {}
            graph[0] = []
            graph[ftn_idx]=[]
            for seq in sequences:
                if len(seq)==0:continue
                if len(seq) == 1:
                    if seq[0] not in graph[0]: # connect the start node with the first node of the sequence
                        graph[0].append(seq[0])
                    if seq[0] not in graph.keys(): # connect the node with target node
                        graph[seq[0]] = [ftn_idx]
                    else:
                        if ftn_idx not in graph[seq[0]]:
                            graph[seq[0]] += [ftn_idx]
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
                        graph[ftn_start] = [ftn_idx]
                    else:
                        if ftn_idx not in graph[ftn_start]:
                            graph[ftn_start] += [ftn_idx]
            return graph

        # A recursive function used by topologicalSort
        def topologicalSortUtil(graph:dict,v,visited:dict, stack):
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
        all_nodes= list(set(all_nodes))+[0,ftn_idx]
        num_nodes=len(all_nodes)

        # build the graph
        graph=get_graph(sequences,ftn_idx)

        # get the path
        # Mark all the vertices as not visited
        visited={}
        for node_idx in all_nodes:
            visited[node_idx]=False

        keys=list(visited.keys())

        stack = []
        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in visited.keys():
            if visited[i] == False:
                topologicalSortUtil(graph,i, visited, stack)
        return stack[1:]  # remove the first element: constructor

    def _get_SV_pairs(self,sv_list:list)->list:
        """
        get SV pairs from a list of SVs
        SV pairs that are made of primitive type SVs have higher priority
        :param sv_list:
        :return:
        """
        index_pairs=[]
        sv_type=[self.FDG.contractInfo.get_sv_type(sv) for sv in sv_list]
        sv_primitive=[]
        sv_non_primitive=[]
        for idx,type in enumerate(sv_type):
            if type.__eq__(fdg.FDG_global.primitive_index):
                sv_primitive.append(idx)
            else:
                sv_non_primitive.append(idx)
        # get SV pairs from primitive types of SVs
        for pair in utils.get_combination_for_a_list(sv_primitive,2):
            index_pairs.append(pair)
        # get SV pairs: one is of primitive type, the other non primitive type
        for pair in utils.get_combination([sv_primitive,sv_non_primitive],2):
            index_pairs.append(pair)
        # get SV pairs from non primitive types of SVs
        for pair in utils.get_combination_for_a_list(sv_non_primitive, 2):
            index_pairs.append(pair)
        sv_pairs=[]
        for pair in index_pairs:
            sv_pairs.append([sv_list[pair[0]],sv_list[pair[1]]])
        return sv_pairs

    def generate_sequences_by_randomly_parentGroups_parents(self, ftn_idx:int)->list:
        """
        1 get parents and group parents based on the SV they write
        2 randomly select parent group subsets,
         randomly select a parent from a parnt group for each parent group subsets
        3 trun parent subset to parent sequence subsets
        merge parent sequence subsets
        :param ftn_idx:
        :param parent_list:
        :return:
        """
        # get parents
        sv_parents = self.FDG.get_parents(ftn_idx)
        if len(sv_parents) == 0: return []
        parent_list = self.FDG.get_parent_list(ftn_idx)
        if len(parent_list)==0: return []

        prt_subset_num_limit=fdg.FDG_global.prt_subset_num_limit

        # randomly select parent subsets
        parent_groups = [value for value in sv_parents.values() if len(value) > 0]

        parent_subsets = []
        if len(parent_groups) == 1:  # consider parent individually
            if len(parent_groups[0]) > prt_subset_num_limit + 1:
                select = np.random.choice(parent_groups[0], size=prt_subset_num_limit, replace=False)
                parent_subsets = [[item] for item in select]
            else:
                parent_subsets = [[item] for item in parent_groups[0]]
        else:
            max_range = 2 ** len(parent_groups)
            if max_range > prt_subset_num_limit + 1:
                select = np.random.choice(list(range(1, max_range)), size=prt_subset_num_limit, replace=False)
            else:
                select = range(1, max_range)
            select_binary = [utils.get_binary(len(parent_groups), value) for value in select]
            for bin_list in select_binary:
                p_subset = []
                for i, bin_ele in enumerate(bin_list):
                    if bin_ele == 1: # means selected
                        # randomly select one parent from the group
                        p_subset.append(np.random.choice(parent_groups[i], size=1, replace=False)[0])
                if p_subset not in parent_subsets:
                    parent_subsets.append(p_subset)


        # trun parent subset to parent sequence subsets and merge if necessary
        final_sequences=[]
        for prt_subset in parent_subsets:

            if len(prt_subset)==1:
                sequences=self.sequenceAndState.get_n_shortest_state_changing_sequences_for_a_function_length(prt_subset[0],1,self.phase1_depth)
                if len(sequences)>0:
                    final_sequences.append(sequences[0]+[ftn_idx])
            else:
                prt_sequence_subset = []
                for prt_idx in prt_subset:
                    # get 1 shortest sequence
                    seq=self.sequenceAndState.get_n_shortest_state_changing_sequences_for_a_function(prt_idx,1)
                    if len(seq)>0:
                        prt_sequence_subset.append(seq[0])
                if len(prt_sequence_subset)>1: # at least should have two sequences to merge
                    sequences=self._get_a_topological_sequence(ftn_idx,prt_sequence_subset)
                    final_sequences.append(sequences)

        return final_sequences

    def generate_sequences_by_randomly_parents(self, ftn_idx:int)->list:
        """
        1 get parents
        2 randomly select parent subsets,

        3 trun parent subset to parent sequence subsets
        merge parent sequence subsets
        :param ftn_idx:
        :param parent_list:
        :return:
        """
        # get parents
        sv_parents = self.FDG.get_parents(ftn_idx)
        if len(sv_parents) == 0: return []
        parent_list = self.FDG.get_parent_list(ftn_idx)
        if len(parent_list)==0: return []

        prt_subset_num_limit=fdg.FDG_global.prt_subset_num_limit


        parent_subsets = []

        max_range = 2 ** len(parent_list)
        if max_range > prt_subset_num_limit + 1:
            select = np.random.choice(list(range(1, max_range)), size=prt_subset_num_limit, replace=False)
        else:
            select = range(1, max_range)
        select_binary = [utils.get_binary(len(parent_list), value) for value in select]
        for bin_list in select_binary:
            p_subset = []
            for i, bin_ele in enumerate(bin_list):
                if bin_ele == 1: # means selected
                    # randomly select one parent from the group
                    p_subset.append(parent_list[i])
            if p_subset not in parent_subsets:
                parent_subsets.append(p_subset)


        # trun parent subset to parent sequence subsets and merge if necessary
        final_sequences=[]
        for prt_subset in parent_subsets:

            if len(prt_subset)==1:
                sequences=self.sequenceAndState.get_n_shortest_state_changing_sequences_for_a_function_length(prt_subset[0],1,self.phase1_depth)
                if len(sequences)>0:
                    final_sequences.append(sequences[0]+[ftn_idx])
            else:
                prt_sequence_subset = []
                for prt_idx in prt_subset:
                    # get 1 shortest sequence
                    seq=self.sequenceAndState.get_n_shortest_state_changing_sequences_for_a_function(prt_idx,1)
                    if len(seq)>0:
                        prt_sequence_subset.append(seq[0])
                if len(prt_sequence_subset)>1: # at least should have two sequences to merge
                    sequences=self._get_a_topological_sequence(ftn_idx,prt_sequence_subset)
                    final_sequences.append(sequences)

        return final_sequences
    def generate_sequences_paper(self,ftn_idx:int)->list:
        """
        get parent combinations
        one combination-> one sequence
        :param ftn_idx:
        :return:
        """
        sv_parents = self.FDG.get_parents(ftn_idx)
        if len(sv_parents) == 0: return []

        sv_list=list(sv_parents.keys())
        # sv_list =[2,3,4]
        generated_sequences=[] # to save the  generated sequences
        sv_combs=[]
        for length in range(1,len(sv_list)+1):
            sv_combs+=utils.get_combination_for_a_list(sv_list,length)
        for sv_comb in sv_combs:
            if len(sv_comb)==1:
                if isinstance(sv_parents, dict):
                    sequences = self.parent_sequences_write_one_SV(ftn_idx, sv_parents[sv_comb[0]],-1,
                                                                   self.phase1_depth)
                    for seq in sequences:
                        if seq not in generated_sequences:
                            generated_sequences.append(seq)
                else:
                    assert False, "data type error"
                continue

            # replace each sv with the corresponding parents
            sv_comb_parents=[sv_parents[sv] for sv in sv_comb]
            parent_combs=utils.get_combination(sv_comb_parents,len(sv_comb))
            for parent_comb in parent_combs:
                parent_sequence_list=[]
                flag_comb=True
                for parent in parent_comb:
                    parent_seq=self.sequenceAndState.get_n_shortest_state_changing_sequences_for_a_function(parent,1)
                    if len(parent_seq)==0:
                        flag_comb=False
                        break # if thre is one parent that does not have a parent sequence, ignore this parent combination
                    parent_sequence_list.append(parent_seq[0])
                if flag_comb:
                    generate_seq=self._get_a_topological_sequence(ftn_idx,parent_sequence_list)
                    if generate_seq not in generated_sequences:
                        generated_sequences.append(generate_seq)

        return generated_sequences

if __name__=="__main__":
    # sequences=[[2,3],[4]]
    sequences = [[4],[2, 3]]
    seqGnt=SequenceGeneration()
    seq=seqGnt._get_a_topological_sequence(5,[[0,1,2],[0,'sv',3,4]])

    print(seq)
