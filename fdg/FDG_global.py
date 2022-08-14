from copy import copy
import numpy as np

global function_coverage_threshold
function_coverage_threshold=100



# max depth of sequence in FDG is set to 5
global phase1_depth_limit
phase1_depth_limit=2


global phase2_include
phase2_include=1  # 1: means include phase 2; others: does not include phase 2;

phase2_method_select=0 #0: p(), 1:p_1(), 2: p_2(),3:p_3()

seq_num_limit=-1 # no limit
seq_len_limit=5


global print_function_coverage
print_function_coverage=0

# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=100

# provide them to FDG_pruner for FDG building
global solidity_path
global contract




# save the coverage (from coverage_plugin)
global coverage
coverage=0

global target_bytecode
target_bytecode=''

# get instruction indices for each function (from soliditycontract)
global ftns_instr_indices
ftns_instr_indices={}

# save the lists that record which instruction is covered (from coverage_plugin)
global ftns_instr_cov
ftns_instr_cov=[]

global mapping
mapping=[]

global solc_indices

global method_identifiers
method_identifiers={}

#===================================
# support executing sequences directly
global sequences
sequences=''

# level: 0: all state variables read
# 1: state variables in conditions
# others: primitive state variables in conditions
global sv_level
sv_level=1
level_phase1=0
level_phase2=1

global primitive_index
primitive_index=2

global non_primitive_index
non_primitive_index=1

