
# support FDG-guided execution and sequence execution
from fdg.FunctionCoverage import FunctionCoverage
from fdg.sequenceExecutionControl import SequenceExecutionControl
from copy import copy, deepcopy
from typing import cast, List
from fdg.node import Node
from fdg import utils
from fdg.FDG import FDG
from fdg.funtion_info import Function_info
from fdg.contractInfopy import ContractInfo
from fdg.sequence import Sequence
from mythril.laser.ethereum.state.world_state import WorldState
from mythril.laser.ethereum.svm import LaserEVM
from mythril.laser.plugin.interface import LaserPlugin
from mythril.laser.plugin.builder import PluginBuilder
from mythril.laser.plugin.plugins import CoveragePluginBuilder
from mythril.laser.plugin.plugins.coverage import coverage_plugin, InstructionCoveragePlugin
from mythril.laser.plugin.plugins.dependency_pruner import get_dependency_annotation, \
     get_ftn_seq_annotation_from_ws
from mythril.laser.plugin.plugins.plugin_annotations import WSDependencyAnnotation, DependencyAnnotation
from mythril.laser.plugin.signals import PluginSkipState
from mythril.laser.ethereum.state.global_state import GlobalState
from mythril.laser.ethereum.transaction.transaction_models import (
    ContractCreationTransaction,
)
from fdg.functionSequence import FunctionSequence

import logging
import fdg.FDG_global
import time
import numpy as np

log = logging.getLogger(__name__)
global contract_ddress
contract_address=0x0
max_number_functions=100


class FDG_prunerBuilder(PluginBuilder):
    name = "fdg-pruner"
    def __call__(self, *args, **kwargs):
        return FDG_pruner(**kwargs)


class FDG_pruner(LaserPlugin):
    """ """
    def __init__(self,instructionCoveragePlugin:InstructionCoveragePlugin):
        """Creates FDG pruner"""
        self._reset()
        self.functionCoverage=FunctionCoverage(instructionCoveragePlugin,\
                                               fdg.FDG_global.ftns_instr_indices,\
                                               fdg.FDG_global.target_bytecode)

    def _reset(self):

        self._iteration_ = 0
        self.solidity = '' # the file path to the solidity file
        self.contract = '' # the contract name to be executed
        self.contract_data=None # save data resulted from Slither

        # coverage related
        self.cov_ftn_instructions_indices = {} # pure names as the keys (limited by the results from compiler)
        self.cov_ftn_identifiers = {} # full names as the keys, directly extracted
        self.cov_ftn_name_to_full_name={} # map pure name to full name
        self.cov_function_coverage={} # use pure or full names as the keys
        self.cov_function_instruction_indices={} # full names as keys
        self.cov_other_instruction_indices={} # pure names as keys ( they are not callable functions)

        # save data during in FDG-guided execution phase
        self.OS_states = {}  # save in-between transaction states

        # save executed sequences
        self.save_no_state_change_sequences = {}
        self.save_state_change_sequences = {}
        self.save_cur_iteration_all_sequences = []  # all sequences in the current iteration
        self.save_cur_iteration_state_change_sequences = []  # all valid sequences in the current iteration


        # related to specifying a set of functions to be executed
        self.fct_hash_2_pc_in_dispatcher={}
        self.instr_list_original=[]
        self.instructions_dict={}

        # phase 2: sequence generation and sequence execution
        self.seq_flag=False
        self.seq_execution=SequenceExecutionControl()
        self.seq_selected_deep_functions=[]
        self.seq_generated_sequences=[]
        self.seq_cur_sequence_index = 0
        self.seq_cur_executing_sequence = []
        self.seq_cur_function_index=0


    def initialize(self, symbolic_vm: LaserEVM) -> None:
        """Initializes the FDG_pruner
        :param symbolic_vm
        """
        self._reset()

        @symbolic_vm.laser_hook("start_sym_exec")
        def start_sym_exec_hook():
            self.solidity = fdg.FDG_global.solidity_path
            self.contract = fdg.FDG_global.contract

            # get contract data
            self.contract_data=ContractInfo(self.solidity, self.contract)
            self.contract_data.get_state_variables_info()
            self.contract_data.get_user_callable_functions_info()

            #
            fdg.FDG_global.ftn_to_idx=self.contract_data.ftn_to_idx
            fdg.FDG_global.idx_to_ftn=self.contract_data.idx_to_ftn
            for key in self.contract_data.functions_info.keys():
                fdg.FDG_global.ftn_to_selector[key]=self.contract_data.functions_info[key]['selector']
                fdg.FDG_global.selector_to_ftn_full_name[self.contract_data.functions_info[key]['selector']]=key



            # get each function's instructions' indices in the whole instruction list
            self.cov_ftn_instructions_indices = fdg.FDG_global.ftns_instr_indices

            for ftn_full_name, identifier in fdg.FDG_global.method_identifiers.items():
                # remve '(...)' from function signature,
                # use pure name as key because self.ftn_instructions_indices uses only pure name as key
                pure_name=str(ftn_full_name).split('(')[0]
                self.cov_ftn_name_to_full_name[pure_name]=ftn_full_name

            for ftn, indices in fdg.FDG_global.ftns_instr_indices.items():
                if ftn in self.cov_ftn_name_to_full_name.keys():
                    self.cov_function_instruction_indices[self.cov_ftn_name_to_full_name[ftn]]=indices
                else:
                    print(f'{ftn} does not have full name')
                    self.cov_other_instruction_indices[ftn] = indices

            # initialize coverage for each function except constructor
            for ftn, ftn_instr_list in fdg.FDG_global.ftns_instr_indices.items():
                # if ftn=='constructor' or ftn=='fallback':continue
                if ftn == 'constructor': continue
                self.cov_function_coverage[self.cov_ftn_name_to_full_name[ftn]]=[0 / len(ftn_instr_list), ftn_instr_list]

        @symbolic_vm.laser_hook("stop_sym_exec")
        def stop_sym_exec_hook():
            if fdg.FDG_global.print_ftn_coverage==1:
                print(f'End of symbolic execution')
                for ftn, ftn_cov in self.cov_function_coverage.items():
                    print("{:.2f}% coverage for {}".format(ftn_cov[0], ftn))



        @symbolic_vm.laser_hook("stop_sym_trans")
        def execute_stop_sym_trans_hook():
            # ----------------------------
            pass

        @symbolic_vm.laser_hook("start_sym_trans_laserEVM")
        def start_sym_trans_hook_laserEVM(laserEVM: LaserEVM):
            """
            ...
            add states to laserEVM.open_states so that they can be used
            as base states in the next iteration of symbolic transaction
            :param laserEVM: instance of LaserEVM
            :return:
            """
            self._iteration_ += 1
            print(f'start: self._iteration_={self._iteration_}')

            if self._iteration_==1:
                # collect the pc for each function hash in dispatcher
                self._collect_pc_for_fct_hashes_in_dispatcher(laserEVM)
                self.save_cur_iteration_all_sequences=[[self.contract_data.idx_to_ftn[idx]] for idx in range(1, len(self.contract_data.idx_to_ftn.keys()))]
                self.save_cur_iteration_state_change_sequences=[]

            # modify instructions
            if self._iteration_ <=fdg.FDG_global.fdg_execution_depth and self._iteration_>1:
                new_states = []
                # modify the instruction list for each open states based on FDG
                for state in laserEVM.open_states:
                    #if not state.constraints.is_possible: continue
                    ftn_seq=get_ftn_seq_annotation_from_ws(state)
                    ftn_name = state.node.function_name

                    children_nodes = self.contract_data.get_children(ftn_name,fdg.FDG_global.sv_level)
                    children_selectors=[node.selector for node in children_nodes]
                    for node in children_nodes:
                        self.save_cur_iteration_all_sequences.append(ftn_seq + [fdg.FDG_global.selector_to_ftn_full_name[node.selector]])

                    # modify the state
                    if len(children_selectors) > 0:
                        modify_state = deepcopy(state)
                        self._modify_dispatcher_in_instruction_list(modify_state, children_selectors)
                        new_states.append(deepcopy(modify_state))

                # update the open states
                laserEVM.open_states = new_states






        @symbolic_vm.laser_hook("stop_sym_trans_laserEVM")
        def stop_sym_trans_hook_laserEVM(laserEVM: LaserEVM):
            """
            - save states
            - only need to states from depth 1 to fdg.FDG_global.depth_all_ftns_reached+1
            - some saved states are used as initial states in sequence execution

            :param laserEVM:
            :return:
            """
            print(f'end: self._iteration_={self._iteration_}')
            if self._iteration_ == 0: return

            # save states
            self._save_states(laserEVM,False)

            # save executed sequences
            print(f'cur_all:{self.save_cur_iteration_all_sequences}')
            print(f'cur_state_change:{self.save_cur_iteration_state_change_sequences}')
            for item in self.save_cur_iteration_state_change_sequences:
                assert len(item)>0
                if item[-1] not in self.save_state_change_sequences.keys():
                    self.save_state_change_sequences[item[-1]]=[item]
                else:
                    if item not in self.save_state_change_sequences[item[-1]]:
                        self.save_state_change_sequences[item[-1]]+=[item]
            for item in self.save_cur_iteration_all_sequences:
                if item not in self.save_cur_iteration_state_change_sequences:
                    if item[-1] not in self.save_no_state_change_sequences.keys():
                        self.save_no_state_change_sequences[item[-1]] = [item]
                    else:
                        if item not in self.save_no_state_change_sequences[item[-1]]:
                            self.save_no_state_change_sequences[item[-1]] += [item]
            # empty saved sequences
            self.save_cur_iteration_state_change_sequences=[]
            self.save_cur_iteration_all_sequences=[]


            # check the code coverage for each function
            if self._iteration_==fdg.FDG_global.fdg_execution_depth:
                self._update_coverage()


            if self.seq_flag:
                if self.seq_execution.function_index==len(self.seq_execution.sequence_cur_in_execution)-1:
                    self._update_coverage()
                self.seq_execution.end_exe_a_function()


            # signal to start sequence execution
            if self._iteration_==fdg.FDG_global.fdg_execution_depth:
                self.seq_flag=True
                laserEVM.open_states=[]


            # sequence execution
            if self.seq_flag:
                # generate the sequences to be executed
                while self.seq_execution.flag_to_generate_sequences:
                    deep_functions = self._get_deep_functions()
                    deep_function_selected = self._select_deep_functions(deep_functions,
                                                                         self.seq_selected_deep_functions)
                    if len(deep_function_selected) > 0:
                        self.seq_selected_deep_functions.append(deep_function_selected)
                        sequences = self._generate_sequences(deep_function_selected)
                        if len(sequences) > 0:
                            self.seq_execution.feed_generated_sequences(sequences)

                    else:
                        # all deep functions are selected once
                        # terminate sequence execution
                        fdg.FDG_global.transaction_count = self._iteration_
                        laserEVM.open_states = []
                        break

                if not self.seq_execution.flag_to_generate_sequences:
                    key, ftn_full_name = self.seq_execution.start_exe_a_function(self.OS_states)
                    if key is not None:
                        laserEVM.open_states = copy(self.OS_states[key])

                    # execute the function
                    if ftn_full_name is not None:
                        new_states = []
                        if ftn_full_name in fdg.FDG_global.ftn_to_selector.keys():
                            ftn_selector = [fdg.FDG_global.ftn_to_selector[ftn_full_name]]
                        else:
                            print(f'function full name: {ftn_full_name}')
                            laserEVM.open_states = []
                            fdg.FDG_global.transaction_count=self._iteration_
                        for state in laserEVM.open_states:
                            # add the sequences to be executed
                            ftn_seq = get_ftn_seq_annotation_from_ws(state)
                            if ftn_seq+[ftn_full_name] not in self.save_cur_iteration_all_sequences:
                                self.save_cur_iteration_all_sequences.append(ftn_seq+[ftn_full_name])

                            modity_state = copy(state)
                            self._modify_dispatcher_in_instruction_list(modity_state, ftn_selector)
                            new_states.append(copy(modity_state))
                        laserEVM.open_states = new_states

        # @symbolic_vm.pre_hook("RETURN")
        # def return_hook(state: GlobalState):
        #     seq = get_dependency_annotation(state).ftn_seq
        #     print(f'return: iteration:{self._iteration_}; ftn:{state.node.function_name}; seq:{seq}')
        #
        # @symbolic_vm.pre_hook("STOP")
        # def stop_hook(state: GlobalState):
        #     seq = get_dependency_annotation(state).ftn_seq
        #     print(f'stop: iteration:{self._iteration_}; ftn:{state.node.function_name}; seq:{seq}')

        # @symbolic_vm.pre_hook("REVERT")
        # def revert_hook(state: GlobalState):
        #     seq = get_dependency_annotation(state).ftn_seq
        #     print(f'revert: iteration:{self._iteration_}; ftn:{state.node.function_name}; seq:{seq}')


        def _transaction_end(state: GlobalState) -> None:
            """
            save function sequences
            :param state:
            """
            # get valid sequences from states
            pass







        @symbolic_vm.laser_hook("add_world_state")
        def world_state_filter_hook(state: GlobalState):
            if isinstance(state.current_transaction, ContractCreationTransaction):
                # Reset iteration variable
                self._iteration_ = 0
                return

    def _print_state_info(state: GlobalState) -> None:
        print(f'==== constraints ====')
        for constraint in state.world_state.constraints:
            print(f'\t {constraint}')
        print(f'==== state.environment.active_account ====')
        print(f'\t {state.environment.active_account.address}')

        print(f'==== storage of the active_account ====')
        for key, value in state.environment.active_account.storage.printable_storage.items():
            print(f'\t key {key}  value {value}')

        print(f'==== memory ====')
        mem_size = state.mstate.memory_size
        for i in range(int(mem_size / 32)):
            word = state.mstate.memory.get_word_at(i)
            print(f'\t {word}')
        print(f'==== stack ====')
        for item in state.mstate.stack:
            print(f'\t {item}')


    def _get_deep_functions(self):
        deep_functions=[]
        for key in self.cov_function_coverage.keys():
            if self.cov_function_coverage[key][0]<=98:
                # only consider callable functions
                if key in self.cov_function_instruction_indices.keys():
                    deep_functions.append(key)
        return deep_functions


    def _select_deep_functions(self,deep_functions:list,selected_deep_functions)->Node:
        ftn_prt_data={}
        for ftn_name in deep_functions:
            if ftn_name in self.contract_data.fdg_parents.keys():
                parents=self.contract_data.fdg_parents[ftn_name]
            else:
                parents=self.contract_data.get_parents(ftn_name,fdg.FDG_global.sv_level)
            for prt in parents:
                if isinstance(prt, Node):
                    if prt.full_name in self.save_state_change_sequences.keys():
                        if ftn_name not in ftn_prt_data.keys():
                            ftn_prt_data[ftn_name]=1
                        else:ftn_prt_data[ftn_name]+=1
        values=list(ftn_prt_data.values())
        if len(values)==0:
            return ""
        if 0 in values:
            values.remove(0)
        min_value=min(values)
        for key,value in ftn_prt_data.items():
            if value==min_value:
                if key not in selected_deep_functions:
                    return key

        return ""


    def _generate_sequences(self,deep_function:str)->list:
        # create the node for the deep_function
        target_ftn_node=Node(fdg.FDG_global.ftn_to_selector[deep_function], \
                             fdg.FDG_global.ftn_to_idx[deep_function], \
                             deep_function,"Not Care")
        if deep_function in self.contract_data.fdg_parents.keys():
            parents=self.contract_data.fdg_parents[deep_function]
        else:
            parents=self.contract_data.get_parents(deep_function,fdg.FDG_global.sv_level)
        ftn_sequence=FunctionSequence(target_ftn_node,parents,\
                                                      self.contract_data.state_variables_info,\
                                                      self.save_state_change_sequences,\
                                                      self.save_no_state_change_sequences[deep_function])

        return ftn_sequence.generate_sequences()


    def _update_coverage(self):
        """
        update coverage and get uncovered functions
        :return:
        """
        instr_cov_record_list = fdg.FDG_global.ftns_instr_cov
        if len(instr_cov_record_list) > 0:
            instr_array = np.array(instr_cov_record_list)
            self.uncovered_functions = []
            self.ftn_special_pc = []
            for ftn, ftn_instr_cov in self.cov_function_coverage.items():
                if ftn_instr_cov[0] == 100: continue
                if ftn in self.cov_function_instruction_indices.keys():
                    status = instr_array[self.cov_function_instruction_indices[ftn]]
                else:
                    status = instr_array[self.cov_other_instruction_indices[ftn]]
                cov_instr = sum(status)
                cov = cov_instr / float(len(status)) * 100
                self.cov_function_coverage[ftn] = [cov, status]

    def _save_states(self,laserEVM:LaserEVM,flag_save_states_1_phase:bool):
        for state in laserEVM.open_states:
            if not state.constraints.is_possible: continue
            ftn_name = state.node.function_name
            ftn_seq=get_ftn_seq_annotation_from_ws(state)
            self.save_cur_iteration_state_change_sequences.append(ftn_seq)


            key=ftn_seq[0]
            for ftn in ftn_seq[1:]:
                key+=";"+ftn
            if key not in self.OS_states.keys():
                self.OS_states[key]=[copy(state)]
            else:
                self.OS_states[key]+=[copy(state)]


    def _collect_pc_for_fct_hashes_in_dispatcher(self, laserEVM:LaserEVM):
        """

        """
        if self._iteration_==1:
            stop_flag=False
            for state in laserEVM.open_states:
                if stop_flag:break # collect only from a state at depth 1
                stop_flag=True
                key = contract_address.value
                code=state.accounts[key].code
                instructions=code.instruction_list
                fct_instr_offsets=[]


                function_hashes=code.func_hashes

                offset_instr=0
                for instruction in instructions:
                    opcode=instruction['opcode']
                    if str(opcode).__eq__('PUSH4'):
                        if instruction['argument'] in function_hashes:
                            if not str(instruction['argument']) in self.fct_hash_2_pc_in_dispatcher.keys():
                                self.fct_hash_2_pc_in_dispatcher[str(instruction['argument'])]=offset_instr
                            else:self.fct_hash_2_pc_in_dispatcher[str(instruction['argument'])]=[offset_instr]+[self.fct_hash_2_pc_in_dispatcher[str(instruction['argument'])]]
                            fct_instr_offsets.append(offset_instr)
                    offset_instr+=1
                    if len(self.fct_hash_2_pc_in_dispatcher)==len(function_hashes):
                        break
                # not yet consider the case of fallback functions
                min_offset=min(fct_instr_offsets)
                max_offst=max(fct_instr_offsets)
                self.instructions_dict['prefix']=instructions[0:min_offset]
                self.instructions_dict['middle'] = instructions[min_offset:max_offst+4]
                self.instructions_dict['suffix'] = instructions[max_offst+4:]

    def _modify_dispatcher_in_instruction_list_old(self,state:WorldState,fct_hashes:list):
        """
            remoeve the code in dispacher that direct the execution flow to functions in fct_hashes
        """

        instr_offsets=[self.fct_hash_2_pc_in_dispatcher[signature] for signature in fct_hashes]
        instr_offsets.sort()
        remove_instr_offsets=[offset for item in instr_offsets for offset in range(item,item+4)]

        max_instr_offset=max(remove_instr_offsets)

        left_instructions=[]
        offset=0
        for instruction in self.instr_list_original[0:max_instr_offset+1]:
            if not offset in remove_instr_offsets:
                left_instructions.append(instruction)
            offset+=1

        left_instructions+=self.instr_list_original[max_instr_offset+1:]

        state.accounts[contract_address.value].code.instruction_list=left_instructions
        function_hashes = state.accounts[contract_address.value].code.func_hashes
        state.accounts[contract_address.value].code.func_hashes=[hash for hash in function_hashes if hash not in fct_hashes]

    def _modify_dispatcher_in_instruction_list(self, state: WorldState, fct_hashes: list):
        """
            remoeve the code in dispacher that direct the execution flow to functions in fct_hashes
        """

        remove_fct_hashes=[ftn_hash for ftn_hash in self.fct_hash_2_pc_in_dispatcher.keys() if ftn_hash not in fct_hashes]
        instr_offsets = [self.fct_hash_2_pc_in_dispatcher[signature] for signature in remove_fct_hashes]

        # Using lambda arguments: expression
        flatten_list = lambda irregular_list:[element for item in irregular_list for element in flatten_list(item)] if type(irregular_list) is list else [irregular_list]
        instr_offsets=flatten_list(instr_offsets)

        instr_offsets.sort()

        remove_instr_offsets = [offset for item in instr_offsets for offset in range(item, item + 5)]

        max_instr_offset = max(remove_instr_offsets)

        left_instructions = []
        offset = len(self.instructions_dict['prefix'])
        for instruction in self.instructions_dict['middle']:
            if not offset in remove_instr_offsets:
                left_instructions.append(instruction)
            else:
                left_instructions.append({"address": instruction["address"], "opcode": "EMPTY"})
            offset += 1

        len_mid=len(left_instructions)
        for i in range(len_mid):
            if left_instructions[len_mid-i-1]['opcode'].__eq__('EMPTY'):
                continue

            if left_instructions[len_mid-i-1]['opcode'].__eq__('DUP1'):
                left_instructions[len_mid-i-1]['opcode']="EMPTY"
            break



        modify_instructions= self.instructions_dict['prefix']+left_instructions+self.instructions_dict['suffix']

        state.accounts[contract_address.value].code.instruction_list = copy(modify_instructions)
        state.accounts[contract_address.value].code.func_hashes = fct_hashes
