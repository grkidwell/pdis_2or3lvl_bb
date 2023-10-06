import numpy as np
import pandas as pd
from winding_temp import dcr_temp

ihlp_filename = 'data/ihlp_core_data.xlsx'

#ckt_params = {'deltaV','duty','fs','Tamb','Idc'}

class Inductor_pdis:
    def __init__(self,ckt_params:dict,Lparams:dict):
        self.ckt = ckt_params
        self.ipp = ckt_params['deltaV']*ckt_params['t_for_deltaV']/Lparams['Lout']/1e-6
        self.idc = ckt_params['Idc']
        self.ind = Lparams
        self.p_core    = self.pcore(ckt_params,Lparams)
        self.tempco = 1/(234.45+25)
        self.t_winding = round(dcr_temp(ckt_params,Lparams,self.p_core,self.tempco),1)
        self.DCR = Lparams['DCR']*(1+self.tempco*(self.t_winding-25))
        self.p_ac = Lparams['K1']*self.ipp**2*ckt_params['fs_phase']**0.5*self.DCR
        self.p_dc = self.DCR*self.idc**2
        self.p_tot = self.p_dc+self.p_ac+self.p_core
        
    def pcore(self,ckt_params:dict,Lparams:dict):
        etckt = ckt_params['deltaV']*ckt_params['t_for_deltaV']/1e-6
        bpk = etckt/Lparams['ET100']*100
        d = ckt_params['duty']; fs = ckt_params['fs_phase']
        fe = fs/(2*np.pi*(d-d**2))
        return Lparams['K0']*fe**(Lparams['Kf']-1)*bpk**Lparams['Kb']*fs*1e-14
    
    def losses(self):
        print(f'Total: {round(self.p_tot,3)}')
        print(f'DC: {round(self.p_dc,3)}')
        print(f'AC: {round(self.p_ac,3)}')
        print(f'core: {round(self.p_core,3)}')
        print(f'Temp: {round(self.t_winding,1)}')


def ind_pdis_obj(ckt_params,lout:float,df_ind_family):
    df_ind = df_ind_family[df_ind_family['Lout']==lout]
    lparams = df_ind.to_dict('records')[0]
    return Inductor_pdis(ckt_params,lparams)

def l_set(ckt_params,family:str,inductor_list:list): #Ldict:dict):
    sheetname = family #Ldict["family"]
    #inductor_list = Ldict["Lout_values"]
    df_ind_family = pd.read_excel(ihlp_filename,sheet_name = sheetname,engine='openpyxl')
    inductor_set = {ind:ind_pdis_obj(ckt_params,ind,df_ind_family) for ind in inductor_list if ind in df_ind_family['Lout'].tolist()}
    pset = [{'Lout':ind_obj.ind['Lout'],
             'Ptot':ind_obj.p_tot,
             'Pdc' :ind_obj.p_dc,
             'Pac' :ind_obj.p_ac,
             'Pcore':ind_obj.p_core,
             'Temp':ind_obj.t_winding
            } for ind_obj in inductor_set.values()
           ]    
    return pd.merge(df_ind_family,pd.DataFrame.from_dict(pset))