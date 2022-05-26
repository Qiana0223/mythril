from slither.core.declarations.function_contract import FunctionContract
from slither.slither import Slither
from fdg.node import Node
import sha3
SV_NOT_CONSIDER=['mapping','array']
class ContractData():
    def __init__(self,solidity_file:str, contract_name:str):
        self.solidity_file=solidity_file
        self.contract_name=contract_name.lstrip().rstrip()

        self.slither_contract=Slither(self.solidity_file).get_contract_from_name(self.contract_name)

        self.state_variables_info = {}
        self.functions_info={}
        self.ftn_to_idx={}
        self.idx_to_ftn={}
        self.fdg_parents={} # save the functions' parent nodes and edge labels
        self.fdg_children={} # save the functions' children nodes and edge labels


    def get_state_variables_info_include_slot(self):
        slot=0
        inner_count=0
        for sv in self.slither_contract.state_variables:
            sv_type=str(sv.type)
            if str(sv_type).startswith('mapping'):
                sv_type='mapping'
            if sv_type in ['string','uint256','mapping','bytes32','address','array']:
                if inner_count>0:
                    slot+=1
                self.state_variables_info[sv.name]=[slot,sv.type]
                print(f'{sv.name}: {slot}: {sv.type}')
                slot+=1
                inner_count = 0
            else:
                if sv_type.__eq__('uint8'):
                    self.state_variables_info[sv.name] = [slot, sv.type]
                    print(f'{sv.name}: {slot}: {sv.type}')
                    inner_count += 1
                else:
                    self.state_variables_info[sv.name] = [None, sv.type]
                    print(f'{sv.name}: {sv.type}: not handled yet')


    def get_state_variables_info(self):
        slot=0
        inner_count=0
        for sv in self.slither_contract.state_variables:
            sv_type=str(sv.type)
            if str(sv_type).startswith('mapping'):
                sv_type='mapping'
                self.state_variables_info[sv.name] = [sv.type,2]
            elif str(sv_type).startswith('array'):
                sv_type = 'array'
                self.state_variables_info[sv.name] = [sv.type,2]
            else:
                self.state_variables_info[sv.name] = [sv.type,1]
            print(f'{sv}:{sv.type}')


    def get_user_callable_functions_info(self):
        ftn_idx=2 # 0: constructor; 1: fallback
        functions_considered=[]

        for f in self.slither_contract.functions:
            if f.name.__eq__('slitherConstructorVariables'): continue
            if f.name.__eq__('slitherConstructorConstantVariables'): continue
            if f.is_constructor: continue
            # only consider public, external functions
            summary = f.get_summary()
            if len(summary) >= 3:
                if summary[2] not in ['public', 'external']:
                    continue

            if f.full_name not in functions_considered:
                functions_considered.append(f.full_name)
                f_info=self._get_a_function_info(f)
                f_info['index']=ftn_idx
                self.functions_info[f.full_name]=f_info
                self.ftn_to_idx[f.full_name] =ftn_idx
                self.idx_to_ftn[ftn_idx] = f.full_name
                ftn_idx+=1

        if 1 not in self.functions_info.keys():
            self.functions_info["fallback"] = {"name": "fallback", "write_sv": [], "read_sv": [], "read_sv_condition": [],
                    "selector": "None","index":1}
            self.ftn_to_idx['fallback']=1
            self.idx_to_ftn[1]='fallback'

        self.functions_info["constructor"]= {"name": "constructor", "write_sv": [], "read_sv": [], "read_sv_condition": [],
                    "selector": "None","index":0}
        self.ftn_to_idx['constructor'] = 0
        self.idx_to_ftn[0] = 'constructor'

    # level: 0: all state variables read
    # 1: state variables in conditions
    # others: primitive state variables in conditions
    def get_children(self, ftn_name:str, level:int)->list:
        children_selector_sv = []
        if ftn_name not in self.functions_info.keys():
            print(f'function with index {ftn_name} does not have static info.')
            return []
        sv_written = self.functions_info[ftn_name]["write_sv"]
        if len(sv_written) == 0: return []
        for sv_w in sv_written:
            for name in self.functions_info.keys()-['constructor']: # do not consider constructor
                sv_read=self._get_sv_read(name,level)
                if sv_w in sv_read:
                    child_node=Node(self.functions_info[name]["selector"],\
                                    self.functions_info[name]["index"],\
                                    name,sv_w)
                    if child_node not in children_selector_sv:
                        children_selector_sv.append(child_node)
        # save children nodes
        if ftn_name not in self.fdg_children.keys():
            self.fdg_children[ftn_name]=children_selector_sv
        return children_selector_sv

    def get_parents(self, ftn_name: str, level: int) -> list:
        parents_selector_sv = []
        if ftn_name not in self.functions_info.keys():
            print(f'function with index {ftn_name} does not have static info.')
            return []
        sv_read=self._get_sv_read(ftn_name,level)
        if len(sv_read) == 0: return []
        for sv_r in sv_read:
            for name in self.functions_info.keys()-['constructor']:
                sv_written = self.functions_info[name]["write_sv"]
                if sv_r in sv_written:
                    prt_node=Node(self.functions_info[name]["selector"], \
                                  self.functions_info[name]["index"],\
                                  name,sv_r)
                    if prt_node not in parents_selector_sv:
                        parents_selector_sv.append(prt_node)
        # save parent nodes
        if ftn_name not in self.fdg_parents.keys():
            self.fdg_parents[ftn_name]=parents_selector_sv
        return parents_selector_sv

    def _get_sv_read(self,ftn_name:str,level:int)->list:
        if level == 1:
            sv_read = self.functions_info[ftn_name]["read_sv_condition"]
        elif level == 0:
            sv_read = self.functions_info[ftn_name]["read_sv"]
        else:
            sv_read = self.functions_info[ftn_name]["read_sv_condition"]
            if len(self.state_variables_info.keys()) > 0:
                sv_read = [sv for sv in sv_read if self.state_variables_info[sv][1] not in SV_NOT_CONSIDER]
        return sv_read

    def _get_a_function_info(self, ftn:FunctionContract):
        w_list = []
        f_w = ftn.all_state_variables_written()
        if len(f_w) > 0:
            w_list = [sv.name for sv in f_w]

        r_list = []
        f_r = ftn.all_state_variables_read()
        if len(f_r) > 0:
            # consider state variables read in all conditions
            r_list = [sv.name for sv in f_r]

        r_list_condition = []
        f_r_condition = ftn.all_conditional_state_variables_read()
        if len(f_r_condition) > 0:
            # consider state variables read in all conditions
            r_list_condition = [sv.name for sv in f_r_condition]

        if ftn.name.__eq__('fallback'):
            func_hash = 'None'
            return {"name":ftn.full_name,"write_sv":w_list,"read_sv":r_list,"read_sv_condition":r_list_condition,"selector":func_hash}
        else:
            func_hash = self._get_function_selector(ftn.full_name)
            return {"name": ftn.full_name, "write_sv": w_list, "read_sv": r_list, "read_sv_condition": r_list_condition,
                    "selector": func_hash}

    def _get_function_selector(self, sig: str) ->str:
        """'
            Return the function id of the given signature
        Args:
            sig (str)
        Return:
            (int)
        """
        s = sha3.keccak_256()
        s.update(sig.encode("utf-8"))
        return '0x'+s.hexdigest()[:8]


if __name__=='__main__':
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/wei_test.sol', 'wei_test')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/wei_test.sol', 'wei_test')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/play_me_quiz.sol', 'play_me_quiz')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/example_contracts/ZetTokenMint.sol', 'ZetTokenMint')
    # ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/example_contracts/AaronTestCoin.sol', 'AaronTestCoin')
    # ftn_info = ContractData('/home/wei/PycharmProjects/Contracts/example_contracts/DxLockEth4Rep.sol', 'Avatar')
    #conData = ContractData('/home/wei/PycharmProjects/Contracts/_wei/Crowdsale.sol', 'Crowdsale')
    #conData = ContractData('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    #conData = ContractData('/media/sf___share_vms/temp/PDC_2.sol', 'PDC_2')

    #conData = ContractData('/media/sf___share_vms/temp/SMT.sol', 'SMT')
    # conData = ContractData('/media/sf___share_vms/temp/Overflow.sol', 'Overflow')
    conData = ContractData('/media/sf___share_vms/temp/PDC_7.sol', 'PDC_7')

    # conData=ContractData('/media/sf___share_vms/temp/PDC.sol', 'PDC')
    conData.get_state_variables_info()
    conData.get_user_callable_functions_info()

    for key, value in conData.functions_info.items():
        print(f'{key}')
        for k,v in value.items():
            print(f'{k}:{v}')
    # ftn_dict= ftn_info.functions_dict_slither()
    # print("===== ftn_dict ====")
    # for key, value in ftn_dict.items():
    #     print("\t{}:  {}".format(key, value))
    # pass



    # a=[24, 29, 189, 34, 118, 194, 265, 39]
    # pairs=get_valid_pc_interval(a,1000)
    # if pc_is_valid(90,pairs):
    #     print(f'90 is in {pairs}')
