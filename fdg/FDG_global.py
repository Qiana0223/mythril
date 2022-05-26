from copy import copy
import numpy as np
global control_level
control_level=2

global seq_num_limit
seq_num_limit=5
global prt_subset_num_limit
prt_subset_num_limit=1

global print_ftn_coverage
print_ftn_coverage=1

# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=2

# provide them to FDG_pruner for FDG building
global solidity_path
global contract

# max depth of sequence in FDG is set to 5
global depth_all_ftns_reached
depth_all_ftns_reached=2

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

#
global ftn_to_idx
ftn_to_idx={}
global idx_to_ftn
idx_to_ftn={}
global ftn_to_selector
ftn_to_selector={}
global selector_to_ftn_full_name
selector_to_ftn_full_name={}

