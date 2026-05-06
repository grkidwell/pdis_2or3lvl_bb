carova_165W = {'vin': 28,
                'vout': 11.1,#13.5,#9, #9,
                'pin': 165,#102.5,#103.7, #164,
                'eff': 0.97,
                'fs':302e3, #367e3, #318e3,   #839k,723k,635k carova at inductor
                'ton_mult':1,
                'tambient':50,
                'tcomponents':{'hs':'tamb', #option for setting to specific temperature.  otherwise will calculate self heating temp at tamb
                               'ls':'tamb',
                               'lout':'tamb'},
                'controller':'rrb48800',#'raa489301', #raa489110, raa489300
                'r_shunt_input':0.01, #10,
                'rboot':5,
                'lout':{'family':'cmle103t',#'hbed053t',#'cmle104t',#'ihlp5050ez01',#'cmll063t', #'hbed053t', #'cmll063t', #'hbed053t', #'ihlp5050ez01', #'hbed053t',
                        'value(uH)':1.0, #2.2, #0.4, #1.0,
                        'config':'single'}, #single, series, parallel
                'lvl_config':'3 level',  #2 level, 4 level
                'hsfet_partnum':'SISH536DN',#'AONR36366', #AONP38324U_HS',#'SIZF5302DT',#'AONS66408', #'SISH536DN', #'RBE030N04', #'SIRA12DDP',#'SISH536DN',#'SIZ342',#'SISH536DN', #'BSZ024N04LS6',#'SISH536DN', #'ISZ0702NLS',#'SIRA74DP',#'SIRA74DP','BSC059N04LS6',#'BSZ063N04LS6',#'SIR426DP',#'SIS488DN',#'SIRA74DP',#'SIS488DN',#'AONS66408',#'AONR66406',#'SIS488DN',#'SISS4410DN',#'SISA14DN',##'SISS4410DN',#'BSC059N04LS6', #'SIS488DN''SIR426DP'
                'lsfet_partnum':'SISS54DN',#'AON7318',#'AONP38324U_LS',#'SIZF5300DT',#'AONS66408', #'SISS54DN',#'BSZ024N04LS6',#'SISA72ADN',#'SISA72ADN',#'AONS66408',#'SISS4410DN',#'SIRA74DP', #'AONS66408',#'BSC059N04LS6',#'BSZ024N04LS6',#'SIRA74DP',
                'q4_partnum':'SISS54DN',#'AON7318', #'AONP36332_HS',#'AON6314', #'SHORT',#'SISS52DN',#'AON6314', #'SHORT', #'BSC020N03MSG',#'BSC020N03MSG',#'SIRA74DP',#,'AON7318',#'AONR36368',#'AON7318',#'SISA14DN',#'SISA14DN', #'SHORT', #'SISS52DN',
                'vgate':10,
                'm_hs':1,
                'm_ls':1,
                #'qrr_vs_i':'sqrt', #'linear', 'constant' = default
                'rd':0.5,
                'caps':{'vin'   :{'partnum':'GRM21BZ71E106KE15L','n':8},
                        'vout'  :{'partnum':'GRM32EC72A106KE05L','n':8},
                        'flying':{'partnum':'C2012X5R1H106K125AC','n':8}},
                'board':'1oz'
               }

jamestown_161W = {'vin': 28,
                'vout': 11.1,#13.5,#9, #9,
                'pin': 161,#102.5,#103.7, #164,
                'eff': 0.97,
                'fs':290e3, #367e3, #318e3,   #839k,723k,635k carova at inductor
                'ton_mult':1,
                'tambient':50,
                'tcomponents':{'hs':'tamb', #option for setting to specific temperature.  otherwise will calculate self heating temp at tamb
                               'ls':'tamb',
                               'lout':'tamb'},
                'controller':'raa489300',#'rrb48800',#'raa489301', #raa489110, raa489300
                'r_shunt_input':0.01, #10,
                'rboot':5,
                'lout':{'family':'cmle103t',#'hbed053t',#'cmle104t',#'ihlp5050ez01',#'cmll063t', #'hbed053t', #'cmll063t', #'hbed053t', #'ihlp5050ez01', #'hbed053t',
                        'value(uH)':1.0, #2.2, #0.4, #1.0,
                        'config':'single'}, #single, series, parallel
                'lvl_config':'3 level',  #2 level, 4 level
                'hsfet_partnum':'SISH536DN',#'AONR36366', #AONP38324U_HS',#'SIZF5302DT',#'AONS66408', #'SISH536DN', #'RBE030N04', #'SIRA12DDP',#'SISH536DN',#'SIZ342',#'SISH536DN', #'BSZ024N04LS6',#'SISH536DN', #'ISZ0702NLS',#'SIRA74DP',#'SIRA74DP','BSC059N04LS6',#'BSZ063N04LS6',#'SIR426DP',#'SIS488DN',#'SIRA74DP',#'SIS488DN',#'AONS66408',#'AONR66406',#'SIS488DN',#'SISS4410DN',#'SISA14DN',##'SISS4410DN',#'BSC059N04LS6', #'SIS488DN''SIR426DP'
                'lsfet_partnum':'SISS54DN',#'AON7318',#'AONP38324U_LS',#'SIZF5300DT',#'AONS66408', #'SISS54DN',#'BSZ024N04LS6',#'SISA72ADN',#'SISA72ADN',#'AONS66408',#'SISS4410DN',#'SIRA74DP', #'AONS66408',#'BSC059N04LS6',#'BSZ024N04LS6',#'SIRA74DP',
                'q4_partnum':'SHORT',#'SISS54DN',#'AON7318', #'AONP36332_HS',#'AON6314', #'SHORT',#'SISS52DN',#'AON6314', #'SHORT', #'BSC020N03MSG',#'BSC020N03MSG',#'SIRA74DP',#,'AON7318',#'AONR36368',#'AON7318',#'SISA14DN',#'SISA14DN', #'SHORT', #'SISS52DN',
                'vgate':10,
                'm_hs':1,
                'm_ls':1,
                #'qrr_vs_i':'sqrt', #'linear', 'constant' = default
                'rd':0.5,
                'caps':{'vin'   :{'partnum':'GRM21BZ71E106KE15L','n':8},
                        'vout'  :{'partnum':'GRM32EC72A106KE05L','n':8},
                        'flying':{'partnum':'C2012X5R1H106K125AC','n':8}},
                'board':'2oz'
               }