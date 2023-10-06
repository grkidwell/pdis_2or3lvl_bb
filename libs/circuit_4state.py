
def up_flag(vin,vout): #duty_under50pct(vin,vout):
    return vout>vin/2

#duty refers to state1/state3 which is 'freewheeling' period / Ts.  For Vo>Vin/2 D=50% results
#in Q1 and Q2 on 100% of time

def duty_freewheel(vin,vout):
    d_up   = (-vin+2*min(vout,vin))/(2*vin)
    d_down = -(-vin+2*vout)/(2*vin)   #not used in actual math.  included for comparison
    return d_up*(-1)**(not up_flag(vin,vout)) #duty_under50pct(vin,vout)


def duty_iramp_up(vin,vout):
    if up_flag(vin,vout):
        return duty_freewheel(vin,vout)
    else:
        return vout/(vin/2)
     

    
def circuit_params(vin,vout,iout,fs,tamb,l_config):
    dv_multiplier   = {'single':1,'series':0.5,'parallel':1}[l_config]
    iout_multiplier = {'single':1,'series':1,'parallel':0.5}[l_config]
    t_state13 = duty_freewheel(vin,vout)/fs
    t_state24 = (1-2*duty_freewheel(vin,vout))/2/fs
    t_Qhs = t_state24+2*t_state13*up_flag(vin,vout)
    t_Qls = t_state24+2*t_state13*(not up_flag(vin,vout))
    eff_est=.95

    return {  'state count':4,
              'vin':vin,
              'vphase':vin/2,
              'deltaV':(vin/2-vout)*dv_multiplier,
              't_for_deltaV':(1-2*duty_freewheel(vin,vout))/2/fs,
              'duty_freewheel':duty_freewheel(vin,vout),
              'duty': duty_iramp_up(vin,vout),
              't_state13': t_state13,
              't_state24': t_state24,
              't_Qhs':t_Qhs,
              't_Qls':t_Qls,
              'fs_phase'  :fs*2,
              'Tamb':tamb,
              'Idc':iout*iout_multiplier,
              'Iinp':iout*iout_multiplier*vout/vin/eff_est
             }
