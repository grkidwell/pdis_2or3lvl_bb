import math
import numpy as np
from numpy import log as ln
from math import e as e
from math import sin as sin
from math import cos as cos
import scipy.optimize as opt
import pandas as pd
import matplotlib.pyplot as plt
import numdifftools as nd



mosfet_filename = r'data/mosfet_data.xlsx'

def get_fet_params(partnumber:str):
    df_fets=pd.read_excel(r'data/mosfet_data.xlsx')
    df = df_fets
    paramdict = dict(df[['parameter',partnumber]].values)
    unitdict = dict(df[['parameter','units']].values)
    scalefactors = dict(df[['parameter','multiplier']].values)
    params = {param:value*scalefactors[param] for param,value in paramdict.items() if unitdict[param] != 'na'}
    params['package']=paramdict['package']
    return params
    
class Fet_cap_vs_vds:
    def __init__(self,fetparams,vds):
        self.fetparams = fetparams
        self.vds    = vds
        self.cgd_0V = 720e-12 #self.fetparams['Ciss_0V']-self.c_gs() #may need to adjust this value if result of function is negative

    def c_gs(self):
        fp=self.fetparams
        return fp['Ciss_Vds2']-fp['Crss_Vds2'] #fp['Ciss_0V']-fp['Crss_1V']*1.2  #20% multiplier for 0V/1V

    def c_gd(self,v_ds:float):
        fp=self.fetparams
        #c_gd_0V = fp['Ciss_0V']-self.c_gs() #may need to adjust this value if result of function is negative
        c_gd_v2 = fp['Crss_Vds2']
        c_gd_1V = fp['Crss_1V']
        a = (1/c_gd_v2-1/self.cgd_0V)
        b = (1/c_gd_1V-1/self.cgd_0V)
        x = math.log(a/b)/math.log(fp['Vds2'])
        c_j2 = 1/(1/c_gd_1V-1/self.cgd_0V)
        return 1/(1/self.cgd_0V+v_ds**x/c_j2)
    
    def c_ds(self,v_ds:float):
        fp=self.fetparams
        c_ds_v2 = fp['Coss_Vds2']-fp['Crss_Vds2']
        c_ds_1V = fp['Coss_1V']-fp['Crss_1V']
        phi = max(1e-19,(fp['Vds2']*c_ds_v2**2-c_ds_1V**2))/(c_ds_1V**2-c_ds_v2**2)
        c_j1 = c_ds_v2*math.sqrt(1+fp['Vds2']/phi)
        return c_j1/math.sqrt(1+v_ds/phi)
    
    def ciss(self,v_ds:float):
        return self.c_gs()+self.c_gd(v_ds)
    
    def coss(self,v_ds:float):
        return self.c_gd(v_ds)+self.c_ds(v_ds)
    
    def crss(self,v_ds:float):
        return self.c_gd(v_ds)
    
    def q_oss(self,vds):  #definite integral
        def f(v_ds:float):
            return self.c_ds(v_ds)+self.c_gd(v_ds)
        f_vectorized = np.vectorize(f)
        v = np.linspace(0,vds,50)
        f_values = f_vectorized(v) 
        return  np.trapz(f_values,v) 
    
    def q_gd(self,vds):  #definite integral
        fp = self.fetparams
        def f(v_ds:float):
            return self.c_gd(v_ds)
        f_vectorized = np.vectorize(f)
        v = np.linspace(fp['Vds_qgd'],vds,50)
        f_values = f_vectorized(v) 
        return  fp['Qgd']+np.trapz(f_values,v) 

class Losses:
    def __init__(self,ckt_params,fs_dcm,*args): 
        self.ic_params,self.hsfet_params,self.lsfet_params,self.vds,self.vgate,idc,ipp,self.m_hs,self.m_ls,rd = args 
        self.lsfp=self.lsfet_params
        self.idc = idc/self.m_ls; self.ipp = ipp/self.m_ls
        self.ckt_params = ckt_params        
        self.state_count = self.ckt_params['state count']
        self.ts = {4:(2*self.ckt_params['t_state13']+2*self.ckt_params['t_state24']),
                   2:self.ckt_params['t_state13']+self.ckt_params['t_state24']}[self.state_count]
        self.fs=fs_dcm*{2:1,4:0.5}[self.state_count]  #1/self.ts    

        self.vth = self.lsfp['Qgs']/self.lsfp['Ciss_Vds2'] 
        #need to tune gate drive to waveform - may need to incorporate datasheet idrivemax
        self.rglsdrvr = {5:self.ic_params['rg_lsdrvr_5V'],
                        10:self.ic_params['rg_lsdrvr_10V']}[self.vgate]
        self.rgls = self.lsfp['Rg']+self.m_ls*self.rglsdrvr        
        self.i_valley = max(0,self.idc-self.ipp/2); self.i_peak = max(self.ipp,self.idc+self.ipp/2)
        
        self.fet_cap = Fet_cap_vs_vds(self.lsfp,self.vds)
        self.summary = {'bd_on':self.bd_f()['on'],
                        'bd_off':self.bd_f()['off'],
                        'cond': self.cond_f(),
                        'ring': self.ring_f(),
                        'gate': self.gate_f()}

    
    def vfwd(self,ifw):
        #need to tune/measure body diode drop vs current
        enabled = False
        vbd_spec = self.lsfet_params['Vbd']
        return vbd_spec+vbd_spec/self.lsfet_params['Id_vbd']*(ifw)**0.5*enabled
    def bd_f(self): #uses 489300 phasenode bd measurements
        v_sgf = 0.1*self.vgate
        t_gsr = self.rgls*self.lsfp['Ciss_0V']*ln(self.vgate/(self.vgate-self.vth))
        t_bd_on = self.ic_params['tsfet_dt_on']#+t_gsr
        t_gsf = self.rgls*self.lsfp['Ciss_0V']*ln(self.vth/v_sgf)
        t_bd_off = t_gsf+self.ic_params['tsfet_dt_off']#+t_gsf+self.rgls*self.lsfp['Ciss_Vds2']*ln(self.vgate/(self.vgate-self.vth))
        return {'on':  self.vfwd(self.i_peak)*self.i_peak*t_bd_on*self.fs,
                'off': self.vfwd(self.i_valley)*self.i_valley*t_bd_off*self.fs,
                'tgsr':t_gsr,
                't_bd_on':t_bd_on,
                'tgsf':t_gsf,
                't_bd_off':t_bd_off} 
                
    def cond_f(self):
        tcoeff = 3500e-6
        tmult = tcoeff*(self.ckt_params['Tamb']-25)
        rdson = {5:self.lsfet_params['Rdson_4.5V'],10:self.lsfet_params['Rdson_10V']}[self.vgate]*(1+tmult)
        t_Qls = self.ckt_params['t_Qls']
        i_fetrms = ((self.idc**2+self.ipp**2/12)*t_Qls/self.ts)**0.5
        return i_fetrms**2*rdson 

    def qrr(self):
            lsp=self.lsfet_params
            qoss = self.fet_cap.q_oss(lsp['Vds_qrr'])
            #old return max((lsp['Qrr']-qoss)/lsp['Id_qrr']*self.i_valley,0)
            #if qoss>qrr then qrr losses aren't counted which makes no sense
            return lsp['Qrr']*(self.i_valley/lsp['Id_qrr'])**0.5

    def ring_f(self):
        qoss_vphase = self.fet_cap.q_oss(self.vds)
        return (self.vds*self.qrr()+qoss_vphase/2*self.vds)*self.fs
        
    def gate_f(self):
        ciss_0V = self.lsfet_params['Ciss_0V']
        qfet_gate = self.vgate*ciss_0V
        vbias = {'no':self.vgate,'yes':self.ckt_params['vin']}[self.ic_params['ldo']]
        return qfet_gate*self.vgate*(vbias/self.vgate)*self.fs
