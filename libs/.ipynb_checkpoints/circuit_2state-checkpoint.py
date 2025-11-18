
def duty(vin,vout):
    return vout/vin
    
def circuit_params(ip):  #only 2 state but using 13 and 24 so keys compatible with circuit_4state.py
    l_config=ip['lout']['config']
    ton_mult = ip['ton_mult']
    vin = ip['vin']; vout=ip['vout']; iout=ip['iout']; fs=ip['fs']/ton_mult; tamb=ip['tambient']
    dv_multiplier   = {'single':1,'series':0.5,'parallel':1}[l_config]
    iout_multiplier = {'single':1,'series':1,'parallel':0.5}[l_config]
    t_state13 = duty(vin,vout)/fs
    t_state24 = (1-duty(vin,vout))/fs
    t_Qhs = t_state13
    t_Qls = t_state24

    eff_est=.95
    
    return {  'state count':2,
              'vin':vin,
              'vphase':vin,
              'deltaV':(vin-vout)*dv_multiplier,
              't_for_deltaV':duty(vin,vout)/fs,
              'duty':duty(vin,vout),
              't_state13': duty(vin,vout)/fs,
              't_state24': (1-duty(vin,vout))/fs,
              't_Qhs':t_Qhs,
              't_Qls':t_Qls,
              'fs_phase'  :fs,
              'ton_mult':ton_mult,
              'Tamb':tamb,
              'Idc':iout*iout_multiplier,
              'Iinp':iout*iout_multiplier*vout/vin/eff_est,
              'ip':ip
             }
