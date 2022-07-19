from mythril.laser.plugin.plugins.coverage import InstructionCoveragePlugin
import numpy as np

class FunctionCoverage():
    def __init__(self,coveragePlugin:InstructionCoveragePlugin,functionInstructionIndices:dict,contractBytecode):
        self.coverage_plugin=coveragePlugin # the plugin that record the instructions visited.
        self.function_instruction_indices= functionInstructionIndices # key: function pure name(only pure name is available from SolidityContract), value: the indices of instructions belonging to the function
        self.contract_bytecode=contractBytecode # the bytecode of the contract that is used to retrieve the instruction coverage status from the coverage plugin

        self.function_coverage={}
        self.contract_coverage=0
        self.ftn_pure_name_to_index={}
        self.deep_function_coverage=[]
        # initialize coverage for each function except constructor
        for ftn, ftn_instr_list in self.function_instruction_indices.items():
            # if ftn=='constructor' or ftn=='fallback':continue
            if ftn == 'constructor': continue
            self.function_coverage[ftn] = 0 / len(ftn_instr_list)

    def set_index_to_ftn_pure_name(self,ftn_full_name_to_index:dict):
        for ftn_full_name,index in ftn_full_name_to_index.items():
            ftn_pure_name=str(ftn_full_name).split('(')[0]
            if ftn_pure_name not in self.ftn_pure_name_to_index.keys():
                self.ftn_pure_name_to_index[ftn_pure_name]=index


    def compute_contract_coverage(self):
        if self.contract_ytecode in self.coverage_plugin.coverage.keys():
            # get the instruction list belonging to the contract
            code_cov = self.coverage_plugin.coverage[self.contract_bytecode ]
            self.coverage = sum(code_cov[1]) / float(code_cov[0]) * 100


    def compute_coverage(self):
        if self.contract_bytecode in self.coverage_plugin.coverage.keys():
            code_cov = self.coverage_plugin.coverage[self.contract_bytecode]
            self.coverage = sum(code_cov[1]) / float(code_cov[0]) * 100 # get contract coverage
            # compute coverage for each function
            instructions_cover_record = code_cov[1]
            if len(instructions_cover_record) > 0:
                instr_array = np.array(instructions_cover_record)
                for ftn,cov in self.function_coverage.items():
                    if cov==100:continue
                    ftn_instruction_indices=self.function_instruction_indices[ftn]
                    status = instr_array[ftn_instruction_indices]
                    cov_instr = sum(status)
                    cov = cov_instr / float(len(status)) * 100
                    self.function_coverage[ftn]=cov


    def get_contract_coverage(self):
        return self.coverage

    def get_function_coverage(self)->dict:
        return self.function_coverage

    def compute_deep_functions(self):
        """
        get a deep functions based the code coverage of each functions

        :return:
        """
        deep_ftn_coverage=[]
        for ftn_pure_name, coverage in self.function_coverage.items():
            if coverage<=98:
                if ftn_pure_name in self.ftn_pure_name_to_index.keys():
                    deep_ftn_coverage.append(self.ftn_pure_name_to_index[ftn_pure_name])
        self.deep_function_coverage=deep_ftn_coverage

    def get_deep_functions_coverage(self): return self.deep_function_coverage
