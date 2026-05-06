# to do - add audible filter as per Lokesh email 3/9/26

# [LK]: The window size is decreased to 3/4th, 2/4th or 1/4th depending on the measured period using an on-chip SM.
# Customers will have two bits control to set the minimum window in their application.
# We don’t plan to go below 1/4th of the window size concerning noise issues. It is decreased to the new window size 
# by counting TBD (12?) clocks for #which period DCM # remains smaller than 50us. If a load transient happens, window 
# will change to larger window size in either 1 or 3 DCM CLK period (trimmable option).

def carova_is_MR(ip) ->bool:
    vin = ip['vin']; vout = ip['vout']; tol = 0.12;
    return (vout>(vin/2*(1-tol))) and (vout<(vin/2*(1+tol)))

def circuit_params(ip):
    ton_mult = ip['ton_mult']
    vin = ip['vin']; vout = ip['vout']; iout=ip['iout']; fs = ip['fs'] ; eff_est= ip['eff']
    tamb=ip['tambient']
    
    lout = ip['lout'];l_config=ip['lout']['config']
    dv_multiplier   = {'single':1,'series':0.5,'parallel':1}[l_config]
    iout_multiplier = {'single':1,'series':1,'parallel':0.5}[l_config]

    if 'ddiv' in ip.keys():
        ddiv = ip['ddiv']
    else:
        ddiv = 2
    d_2state = vout/vin
    ts=1/fs/2  
    tstate = {1:d_2state/ddiv*ts}
    tstate[2] = tstate[1]*2*(ddiv-1)  
    tstate[3] = ts - tstate[1]-tstate[2]
    tstate = {state:round(t,10) for state,t in tstate.items()}

    deltaV = {1:vin-vout, 2:vin/2-vout, 3:-vout}
    deltaVperinductor = {state:dv*dv_multiplier for state,dv in deltaV.items()}
    voltsecmax = max([abs(deltaVperinductor[state]*tstate[state]) for state in tstate.keys()])
    t_Qhs = 2*tstate[1]+tstate[2]
    t_Qls = tstate[2]+2*tstate[3]
    
    return {  'state count':6,
              'vin':vin,
              'vphase':vin/2,  #used to determine fet losses
              'deltaV':deltaVperinductor,  #dictionary for 3states
              't_state':tstate, #dictionary for 3states
              'volt-sec':voltsecmax,  #per inductor
              't_Qhs':t_Qhs,
              't_Qls':t_Qls,
              'fs_phase'  :fs*2,
              'ton_mult':ton_mult,
              'Tamb':tamb,
              'Idc':iout*iout_multiplier,
              'Iinp':iout*iout_multiplier*vout/vin/eff_est,
              'ip':ip
             }
