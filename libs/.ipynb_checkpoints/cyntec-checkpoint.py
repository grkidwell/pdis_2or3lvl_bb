#in this version, will remove lparams dependent variable from Inductor_pdiss constructor 
#and create lparams dictionary within same constructor
#also change lset function 

import numpy as np
import pandas as pd
from winding_temp_dcm import dcr_temp

from circuit_6state import circuit_params as circuit_params_6state
from circuit_4state import circuit_params as circuit_params_4state
from circuit_2state import circuit_params as circuit_params_2state

cyntec_filename = r'data/cyntec_inductor_data.xlsx'

class Inductor_pdis:
    def __init__(self,inp_params): 
        self.ip = inp_params.copy()
        self.ts = 1/self.ip['fs']
                
        self.update_ckt_params()
        
        l_ip = self.ckt['ip']['lout']
        self.ind = lparams(l_ip['value(uH)'],create_ind_family_df(l_ip['family']))   #Lparams.copy()
        self.idc = self.ckt['Idc']
                
        self.do_dcm_stuff()
        
        self.ton_mult=self.ckt['ton_mult']
        
        #phasenode dcm and ccm frequency
        #self.fs_dcm = round(self.dcm_ratio/(self.ts/2)) 
        
        self.ind['K1']=0  #no ac winding loss.  see winding_temp.py
        self.p_core    = self.pcore()
        
        self.tempco = 1/(234.45+25)
        if ('tcomponents' in self.ip and
            'lout' in self.ip['tcomponents'] and
            isinstance(self.ip['tcomponents']['lout'],int)):
            self.t_winding = self.ip['tcomponents']['lout']
        else:
            self.t_winding = round(dcr_temp(self.irms_dcm,self.ckt,self.ind,self.p_core,self.tempco),1)
        self.DCR = self.ind['DCR']*(1+self.tempco*(self.t_winding-25))
        self.p_dc = self.DCR*self.irms_dcm**2
        
        self.p_tot = self.p_dc+self.p_core
        self.summary = {'dcr':self.p_dc,
                        'core':self.p_core,
                        'ipp':self.ipp,
                        'fs_dcm':self.fs_dcm,
                        'ton_mult':self.ton_mult,
                        'irms_dcm':self.irms_dcm}

    def update_ckt_params(self):  #can re-rerun this function after change self.ip['lvl_config']
        def carova_is_MR() ->bool:
            vin = self.ip['vin']; vout = self.ip['vout']; tol = 0.12;
            return (vout>(vin/2*(1-tol))) and (vout<(vin/2*(1+tol)))

        if 'carova' not in self.ip['lvl_config']: #lvl_config options - '2 level', '3 level', '3 level carova'
            self.ckt = {'2 level':circuit_params_2state(self.ip),
                        '3 level':circuit_params_4state(self.ip)}[self.ip['lvl_config']]
        else:
            self.ckt = {True: circuit_params_6state(self.ip),
                        False:circuit_params_4state(self.ip)}[carova_is_MR()]


    def irms_fsdcm_carova_MR(self):
        tstate = self.ckt['t_state'];deltaV = self.ckt['deltaV']
        ipp = {state:round(dv*tstate[state]/self.ind['Lout']/1e-6,3) for state,dv in deltaV.items()}   
        ipulse_avg_per_segment = {1:ipp[1]/2,2:ipp[1]+ipp[2]/2,3:abs(ipp[3])/2}
        ipulse_avg = sum([time*ipulse_avg_per_segment[state] for state,time in tstate.items()])/sum(tstate.values())
        dcm_ratio = min(1,self.idc/ipulse_avg )
        fs_dcm = dcm_ratio/(self.ts/2)
        ipulse_start = {1:0, 2:ipp[1], 3:ipp[1]+ipp[2]}
        i_lout_points = {0:self.idc - ipulse_avg,
                         1:ipulse_start[2]-ipulse_avg+self.idc,
                         2:ipulse_start[3]-ipulse_avg+self.idc,
                         3:self.idc - ipulse_avg}
        ilp = i_lout_points
        def i_ms_segment(idx1,idx2):
            return (ilp[idx1]**2+ilp[idx1]*ilp[idx2]+ilp[idx2]**2)/3  #from Erickson p.748
        self.i_start_bystate = {time+1:current for time,current in ilp.items() if time<3}
        self.i_stop_bystate = {time:current for time,current in ilp.items() if time>0}
        self.i_ms_bystate = {state:i_ms_segment(state-1,state) for state,time in tstate.items()}
        i_lout_rms = (sum([time*self.i_ms_bystate[state] for state,time in tstate.items()])*fs_dcm)**0.5        
        #i_lout_rms = (sum([time*i_ms_segment(state-1,state) for state,time in tstate.items()])*fs_dcm)**0.5
        return i_lout_rms, fs_dcm

    def irms_fsdcm_2or4state(self):
        tstate = self.ckt['t_state']
        is2state = self.ckt['state count']==2
        ipp = self.ipp #self.ckt['volt-sec']/self.ind['Lout']/1e-6
        ipulse_avg=ipp/2
        dcm_ratio = min(1,self.idc/ipulse_avg )
        fs_dcm = round(dcm_ratio/(self.ckt['t_state13']+self.ckt['t_state24']),0) #dcm_ratio/(self.ts/2)
        d1 = self.ckt['t_state13']*fs_dcm
        d2 = self.ckt['t_state24']*fs_dcm
        self.i_start_bystate = {1:self.idc+1/2*ipp*(-1)**(is2state or (self.ckt['d_up_flag'])),
                                2:self.idc+1/2*ipp*(-1)**(not(is2state or self.ckt['d_up_flag']))}
        def i_ms_segment(i1,i2):
            return (i1**2+i1*i2+i2**2)/3  #from Erickson p.748   
        self.i_ms_bystate = {1:i_ms_segment(self.i_start_bystate[1],self.i_start_bystate[2]),
                             2:i_ms_segment(self.i_start_bystate[2],self.i_start_bystate[1])}
        self.i_lout_rms_check =(sum([time*self.i_ms_bystate[state] for state,time in tstate.items()])*fs_dcm)**0.5 
        i_lout_rms = (self.idc**2+(d1+d2)/12*self.ipp**2)**0.5
        return i_lout_rms, fs_dcm
        
    def do_dcm_stuff(self):
        self.ipp = self.ckt['volt-sec']/self.ind['Lout']/1e-6
        if self.ckt['state count'] == 6:
            self.irms_dcm, fs_dcm = self.irms_fsdcm_carova_MR()
        else:
            self.irms_dcm, fs_dcm = self.irms_fsdcm_2or4state()
            
        ccm_hyst  = 1.1
        if fs_dcm > 2*self.ip['fs']:
            ip=self.ckt['ip']; ip['ton_mult']= 1 #this also updates self.ckt['ip']['ton_mult'] to 1
            self.update_ckt_params()
            self.ipp = self.ckt['volt-sec']/self.ind['Lout']/1e-6
            dcm_ratio=1            
        #self.ipp = ipp
        #self.dcm_ratio = dcm_ratio #min(1,dcm_ratio)
        self.fs_dcm = fs_dcm
            
    def pcore(self): 
        return self.ind['Ka']*self.fs_dcm**(self.ind['Kx'])*(self.ind['Kb']*self.ipp)**self.ind['Ky']
    

    
    def losses(self):
        print(f'Total: {round(self.p_tot,3)}')
        print(f'DC: {round(self.p_dc,3)}')
        #print(f'AC: {round(self.p_ac,3)}')
        print(f'core: {round(self.p_core,3)}')
        print(f'Temp: {round(self.t_winding,1)}')

def create_ind_family_df(familyname:str):
    return pd.read_excel(cyntec_filename,sheet_name = familyname,engine='openpyxl')
    
def closest_value(input_list, input_value):
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return arr[i]

def lparams(lout:float,df_ind_family):
    lout_column = df_ind_family['Lout']
    cv = closest_value(lout_column.tolist(),lout)
    df_ind = df_ind_family[lout_column==cv]
    return df_ind.to_dict('records')[0]

#this next function should be unnecessary and depricated
def ind_pdis_obj(ckt_params,lout:float,df_ind_family):
    cp=ckt_params.copy()
    lout_column = df_ind_family['Lout']
    cv = closest_value(lout_column.tolist(),lout)
    df_ind = df_ind_family[lout_column==cv]
    lparams = df_ind.to_dict('records')[0]
    return Inductor_pdis(cp,lparams)

#def l_set(ckt_params,ind_family:str,inductor_list:list): #Ldict:dict):
def l_set(inp_params,ind_family:str,inductor_list:list): #Ldict:dict):
    ip = inp_params.copy()
    ip['lout']['family'] = ind_family
    # cp = {'2 level':circuit_params_2state(ip),
    #           '3 level':circuit_params_4state(ip)}\
    #           [ip['lvl_config']]

    #cp=ckt_params.copy()
    #cp['ip']['lout']['family'] = ind_family
    df_ind_family = create_ind_family_df(ind_family)
    #sheetname = family #Ldict["family"]
    #inductor_list = Ldict["Lout_values"]
    #df_ind_family = create_ind_family_df(family) #pd.read_excel(cyntec_filename,sheet_name = sheetname,engine='openpyxl')
    def ip_replace_ind(ind):
        ip['lout']['value(uH)'] = ind
        return ip
    #inductor_set = {ind:ind_pdis_obj(cp,ind,df_ind_family) for ind in inductor_list if ind in df_ind_family['Lout'].tolist()}
    inductor_set = {ind:Inductor_pdis(ip_replace_ind(ind)) for ind in inductor_list if ind in df_ind_family['Lout'].tolist()}
    pset = [{'Lout':ind_obj.ind['Lout'],
             'Ptot':ind_obj.p_tot,
             'Pdc' :ind_obj.p_dc,
             'Pcore':ind_obj.p_core,
             'Temp':ind_obj.t_winding,
             'Fs':ind_obj.fs_dcm
            } for ind_obj in inductor_set.values()
           ]    
    return pd.merge(df_ind_family,pd.DataFrame.from_dict(pset))