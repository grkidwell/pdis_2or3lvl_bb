# 3lvl carova fets in vmid mode (6state) may transition at multiple states
# look up the corresponding switch current in the i_start_by_state dictionary 
# which is attached to the inductor object

states = {'2' : [1,2],
          '4L': [1,2],
          '4H': [1,2],
          '6' : [1,2,3]}

flycap_states =  {'2' : [0],
                  '4L': [2],
                  '4H': [2],
                  '6' : [2]}

def flycap_cond_states(cp):
    lvl_config = cp['ip']['lvl_config']
    state_count = str(cp['state count'])
    if state_count == '4':
        state_count = state_count+{False:'L',True:'H'}[cp['d_up_flag']]
    return flycap_states[state_count]