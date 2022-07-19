from copy import copy
import numpy as np
global control_level
control_level=2

global seq_num_limit
seq_num_limit=5

global num_sequences_write_1_SV
num_sequences_write_1_SV=2


global prt_subset_num_limit
prt_subset_num_limit=5

global print_ftn_coverage
print_ftn_coverage=1

# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=50

# provide them to FDG_pruner for FDG building
global solidity_path
global contract

# max depth of sequence in FDG is set to 5
global fdg_execution_depth
fdg_execution_depth=2

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
