
# dictionary values are the states at transition.  
# 3lvl carova fets in vmid mode (6state) may transition at multiple states
# look up the corresponding switch current in the i_start_by_state dictionary 
# which is attached to the inductor object

states = {'2' : [1,2],
          '4L': [1,2],
          '4H': [1,2],
          '6' : [1,2,3]}

#I don't know how useful below function is, but wrote it anyway
#if used, it will need to be fixed because it generates [1] for complement of [2,3] instead of [1,1].
def state_complement(fetstatedict):
    return {key:[element for element in states[key] if element not in value] for key,value in fetstatedict.items()}

#these are actual states for each of the topologies/modes, which is all you need
#switch FETs
Q3Q4on = {'2' : [1],
          '4L': [2],
          '4H': [1],
          '6' : ['1a','1b']}

Q3Q4off = {'2' : [2],
           '4L': [1],
           '4H': [2],
           '6' : [2,3]}

Q3Q4cond = {'2' : [1],
            '4L': [2],
            '4H': ['1a',2,'1b'],
            '6' : ['1a',2,'1b']}

#SR FETs
Q1Q2on =  {'2' : [2],
           '4L': [1],
           '4H': [2],
           '6' : [2,3]}

Q1Q2off = {'2' : [1],
           '4L': [2],
           '4H': [1],
           '6' : ['1a','1b']}

Q1Q2cond = {'2' : [2],
            '4L': ['1a',2,'1b'],
            '4H': [2],
            '6' : [2,'3a','3b']}

def hs_switch_states(cp):
    lvl_config = cp['ip']['lvl_config']
    state_count = str(cp['state count'])
    if state_count == '4':
        state_count = state_count+{False:'L',True:'H'}[cp['d_up_flag']]
    return {'on':Q3Q4on[state_count],
            'off':Q3Q4off[state_count]}

def hs_cond_states(cp):
    lvl_config = cp['ip']['lvl_config']
    state_count = str(cp['state count'])
    if state_count == '4':
        state_count = state_count+{False:'L',True:'H'}[cp['d_up_flag']]
    return Q3Q4cond[state_count]

def ls_switch_states(cp):
    lvl_config = cp['ip']['lvl_config']
    state_count = str(cp['state count'])
    if state_count == '4':
        state_count = state_count+{False:'L',True:'H'}[cp['d_up_flag']]
    return {'on':Q1Q2on[state_count],
            'off':Q1Q2off[state_count]}

def ls_cond_states(cp):
    lvl_config = cp['ip']['lvl_config']
    state_count = str(cp['state count'])
    if state_count == '4':
        state_count = state_count+{False:'L',True:'H'}[cp['d_up_flag']]
    return Q1Q2cond[state_count]
