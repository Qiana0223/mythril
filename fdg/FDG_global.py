from copy import copy
import numpy as np

global function_coverage_threshold
function_coverage_threshold=100



# max depth of sequence in FDG is set to 5
global phase1_depth_limit
phase1_depth_limit=2


global flag_phase2
flag_phase2=1  # 1: means include phase 2; others: does not include phase 2;



seq_num_limit=0



global print_ftn_coverage
print_ftn_coverage=1

# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=50

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


global primitive_index
primitive_index=2

global non_primitive_index
non_primitive_index=1
