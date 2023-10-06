import numpy as np
import pandas as pd
from winding_temp_dcm import dcr_temp

cyntec_filename = r'data/cyntec_inductor_data.xlsx'

#ckt_params = {'deltaV','duty','fs','Tamb','Idc'}

class Inductor_pdis:
    def __init__(self,ckt_params:dict,Lparams:dict):
        self.ckt = ckt_params.copy()
        self.ipp = self.ckt['deltaV']*self.ckt['t_for_deltaV']/Lparams['Lout']/1e-6
        self.idc = self.ckt['Idc']
        
        self.dcm_ratio=min(1,self.idc/(self.ipp/2))
        
        #phasenode dcm and ccm frequency
        self.fs_dcm = self.dcm_ratio/(self.ckt['t_state13']+self.ckt['t_state24'])
        
        self.irms_dcm = self.i_rms_dcm()
        
        
        self.ind = Lparams.copy()
        self.ind['K1']=0  #no ac winding loss.  see winding_temp.py
        self.p_core    = self.pcore()#self.ckt,Lparams)
        self.tempco = 1/(234.45+25)
        self.t_winding = round(dcr_temp(self.irms_dcm,self.ckt,self.ind,self.p_core,self.tempco),1)
        self.DCR = self.ind['DCR']*(1+self.tempco*(self.t_winding-25))
        #self.p_ac = Lparams['K1']*self.ipp**2*ckt_params['fs']**0.5*self.DCR
        #self.p_dc = self.DCR*self.idc**2
        self.p_dc = self.DCR*self.irms_dcm**2
        self.p_tot = self.p_dc+self.p_core
        self.summary = {'dcr':self.p_dc,
                        'core':self.p_core}
        
    def pcore(self): 
        fs = self.fs_dcm
        return self.ind['Ka']*fs**(self.ind['Kx'])*(self.ind['Kb']*self.ipp)**self.ind['Ky']
    
    def i_rms_dcm(self):
        d1 = self.ckt['t_state13']*self.fs_dcm
        d2 = self.ckt['t_state24']*self.fs_dcm
        return (self.idc**2+(d1+d2)/12*self.ipp**2)**0.5
    
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

def ind_pdis_obj(ckt_params,lout:float,df_ind_family):
    cp=ckt_params.copy()
    lout_column = df_ind_family['Lout']
    cv = closest_value(lout_column.tolist(),lout)
    df_ind = df_ind_family[lout_column==cv]
    lparams = df_ind.to_dict('records')[0]
    return Inductor_pdis(cp,lparams)

def l_set(ckt_params,df_ind_family,inductor_list:list): #Ldict:dict):
    cp=ckt_params.copy()
    #sheetname = family #Ldict["family"]
    #inductor_list = Ldict["Lout_values"]
    #df_ind_family = create_ind_family_df(family) #pd.read_excel(cyntec_filename,sheet_name = sheetname,engine='openpyxl')
    inductor_set = {ind:ind_pdis_obj(cp,ind,df_ind_family) for ind in inductor_list if ind in df_ind_family['Lout'].tolist()}
    pset = [{'Lout':ind_obj.ind['Lout'],
             'Ptot':ind_obj.p_tot,
             'Pdc' :ind_obj.p_dc,
             'Pcore':ind_obj.p_core,
             'Temp':ind_obj.t_winding
            } for ind_obj in inductor_set.values()
           ]    
    return pd.merge(df_ind_family,pd.DataFrame.from_dict(pset))