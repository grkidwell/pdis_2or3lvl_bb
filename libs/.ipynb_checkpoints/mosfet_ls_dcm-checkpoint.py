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

from mosfet_iswitch_states import ls_switch_states, ls_cond_states 




mosfet_filename = r'data/mosfet_data.xlsx'

def get_fet_params(partnumber:str):
    #df_fets=pd.read_excel(r'data/mosfet_data.xlsx')
    df_fets=pd.read_excel(r'data/mosfet_data.xlsx', sheet_name='transpose')
    #yes, we are tranposing the transposed sheet.  didn't feel like changing the rest of the code
    df3 = df_fets.transpose().reset_index()#drop=True)       df2 doesn't exist and is legacy from test function
    column_names = df3.iloc[0].values.tolist()
    df3.columns=column_names
    df3.drop(index=df3.index[0], axis=0, inplace=True)
    df3.reset_index(inplace=True)
    df3.drop(['index'],axis=1,inplace=True)
    #df3[['multiplier']]=df3[['multiplier']].applymap('{:.0E}'.format)
    df = df3 #df_fets
    paramdict = dict(df[['parameter',partnumber]].values)
    unitdict = dict(df[['parameter','units']].values)
    scalefactors = dict(df[['parameter','multiplier']].values)
    params = {param:value*scalefactors[param] for param,value in paramdict.items() if unitdict[param] != 'na'}
    params['package']=paramdict['package']
    return params

def total_currents(lout_obj):
    ip = lout_obj.ckt['ip']
    idc = lout_obj.idc*{'single':1,'series':1,'parallel':2}[ip['lout']['config']]
    ipp = lout_obj.ipp*{'single':1,'series':1,'parallel':2}[ip['lout']['config']]
    return {'idc':idc,'ipp':ipp}

def stripint(a):
    if type(a) == str:
        a = int(a[0])
    return a

class Fet_cap_vs_vds:
    def __init__(self,fetparams,vds):
        self.fetparams = fetparams
        self.vds    = vds
        self.cgd_0V = max(self.fetparams['Crss_1V']+100e-12,self.fetparams['Ciss_0V']-self.c_gs()) #720e-12 #may need to adjust this value if result of function is negative

    def c_gs(self):
        fp=self.fetparams
        return  max(220e-12,fp['Ciss_0V']-fp['Crss_1V'])  #20% multiplier for 0V/1V,   or fp['Ciss_Vds2']-fp['Crss_Vds2']

    def c_gd(self,v_ds:float):
        fp=self.fetparams
        #c_gd_0V = fp['Ciss_0V']-self.c_gs() #may need to adjust this value if result of function is negative
        c_gd_v2 = fp['Crss_Vds2']
        c_gd_1V = fp['Crss_1V']
        a = (1/c_gd_v2-1/self.cgd_0V)
        
        b = max(1e8,(1/c_gd_1V-1/self.cgd_0V))
        
        ratio = a/b #max(2,a/b)
        x = math.log(ratio)/math.log(fp['Vds2'])
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
    def __init__(self,hs_losses_obj): 
        self.lout_obj = hs_losses_obj.lout_obj;self.ckt_params = hs_losses_obj.ckt_params        
        self.ip = self.ckt_params['ip']
        
        self.vds=self.ckt_params['vphase'];self.vgate=self.ip['vgate']
        self.m_hs=self.ip['m_hs'];self.m_ls=self.ip['m_ls'];self.rd=self.ip['rd']
        self.ic_params = hs_losses_obj.ic_params
        self.hsfet_params = hs_losses_obj.hsfet_params
        self.lsfet_params = hs_losses_obj.lsfet_params       
        
        self.state_count = self.ckt_params['state count']
        self.fs=self.lout_obj.fs_dcm*{2:1,4:0.5,6:0.5}[self.state_count]

        self.vth = self.lsfet_params['Qgs']/self.lsfet_params['Ciss_Vds2'] 

        self.fet_cap = Fet_cap_vs_vds(self.lsfet_params,self.vds)
        
        self.ls_sw_states = ls_switch_states(self.ckt_params)
        self.i_scaler = {'single':1,'series':1,'parallel':2}[self.ip['lout']['config']]/self.m_ls        

        self.bd_on_loss_dict = {state:self.bd_f('on',state) for state in self.ls_sw_states['on']}
        self.bd_off_loss_dict = {state:self.bd_f('off',state) for state in self.ls_sw_states['off']}
        self.ring_loss_dict = {state:self.ring_f(state) for state in self.ls_sw_states['on']}

        self.ls_cond_states = ls_cond_states(self.ckt_params)
        self.fetrms_dcm = self.rms_dcm_calculate()
        
        self.summary = {'bd_on':sum(self.bd_on_loss_dict.values()),
                        'bd_off':sum(self.bd_off_loss_dict.values()),
                        'ring': sum(self.ring_loss_dict.values()),
                        'gate': self.gate_f()}
        if ('tcomponents' in self.ip and
            'ls' in self.ip['tcomponents'] and
            isinstance(self.ip['tcomponents']['ls'],int)):
            self.temp = self.ip['tcomponents']['ls']
        else:
            self.temp = self.temp_f()
        self.summary['cond'] =  self.cond_f()
                        

    

    
        
    def bd_f(self,on_or_off,state): #uses 489300 phasenode bd measurements
        state = stripint(state)
        def vfwd(ifw):
            #need to tune/measure body diode drop vs current
            enabled = False
            vbd_spec = self.lsfet_params['Vbd']
            return vbd_spec+vbd_spec/self.lsfet_params['Id_vbd']*(ifw)**0.5*enabled
        def pfwd(itot,time):
            ifw = itot*self.i_scaler  #for fets in parallel
            return vfwd(ifw)*ifw*time*self.fs

        rglsdrvr = {5:self.ic_params['rg_lsdrvr_5V'],
                        10:self.ic_params['rg_lsdrvr_10V']}[self.vgate]
        rgls = self.lsfet_params['Rg']+self.m_ls*rglsdrvr 
        v_sgf = 0.1*self.vgate
        t_gsr = rgls*self.lsfet_params['Ciss_0V']*ln(self.vgate/(self.vgate-self.vth))
        t_bd_on = self.ic_params['tsfet_dt_on']#+t_gsr
        t_gsf = rgls*self.lsfet_params['Ciss_0V']*ln(self.vth/v_sgf)
        t_bd_off = t_gsf+self.ic_params['tsfet_dt_off']#+t_gsf+self.rgls*self.lsfet_params['Ciss_Vds2']*ln(self.vgate/(self.vgate-self.vth))
        t_bd = {'on':t_bd_on,
                'off':t_bd_off}[on_or_off]
        itot = self.lout_obj.i_start_bystate[state]*self.i_scaler
        return pfwd(itot,t_bd)
        
        # return {'on':  sum(pbd_on_dict.values()),
        #         'off': sum(pbd_off_dict.values()),
        #         'tgsr':t_gsr,
        #         't_bd_on':t_bd_on,
        #         'tgsf':t_gsf,
        #         't_bd_off':t_bd_off} 
                
    def rms_dcm_calculate(self):
        ts = 1/self.fs
        i_scaler = {'single':1,'series':1,'parallel':2}[self.ip['lout']['config']]/self.m_ls
        #period_multiplier = {2:1,4:2,6:2}[self.state_count]
        cond_states = self.ls_cond_states
        ims = self.lout_obj.i_ms_bystate
        tstate = self.ckt_params['t_state']
        t_fullcycle = 1/self.fs
        ms_tot = sum([tstate[stripint(state)]*ims[stripint(state)] for state in cond_states])
        return (ms_tot/t_fullcycle)**0.5*i_scaler

    def temp_f(self):
        pfixed = sum([val for key,val in self.summary.items() if key in ['bd_on','bd_off','ring']])
        Rth = self.lsfet_params['RthJA']
        tamb = self.ckt_params['Tamb']
        tempco = 3500e-6
        rdson_25C = {5:self.lsfet_params['Rdson_4.5V'],10:self.lsfet_params['Rdson_10V']}[self.vgate]
        rdson_term = self.fetrms_dcm**2*rdson_25C
        return abs(((25*tempco-1)*rdson_term-tamb-pfixed*Rth)/(rdson_term*tempco-1))
        
    def cond_f(self):
        tcoeff = 3500e-6
        tmult = tcoeff*(self.temp-25)
        rdson = {5:self.lsfet_params['Rdson_4.5V'],10:self.lsfet_params['Rdson_10V']}[self.vgate]*(1+tmult)
        return self.fetrms_dcm**2*rdson 

    def qrr(self,state,*args:str):
        state = stripint(state)
        i_bd = self.lout_obj.i_start_bystate[state]*self.i_scaler
        lsp=self.lsfet_params
        qrr_ls = lsp['Qrr']
        qoss = self.fet_cap.q_oss(lsp['Vds_qrr'])
        if 'print' in args:
            print(f'qrr: {qrr_ls}')
            print(f'qoss: {qoss}')
        if qrr_ls > qoss:
            qrr_net = qrr_ls-qoss
        else:
            qrr_net = qrr_ls
        qrr_net = qrr_ls
        if 'qrr_vs_i' not in self.ip.keys():
            self.ip['qrr_vs_i'] = 'constant'
        qrr_exp = {'constant':0,'linear':1,'sqrt':0.5}[self.ip['qrr_vs_i']]
        qrr_net = qrr_net*(i_bd/lsp['Id_qrr'])**qrr_exp
        # if 'linear' in args:
        #     qrr_net = qrr_net*self.i_valley/lsp['Id_qrr']
        # elif 'sqrt' in args:
        #     qrr_net = qrr_net*(self.i_valley/lsp['Id_qrr'])**0.5
        return max(qrr_net,0)
        #if qoss>qrr then qrr losses aren't counted which makes no sense

    def ring_f(self,state):   
        qoss_vphase = self.fet_cap.q_oss(self.vds)
        return (self.vds*self.qrr(state)+qoss_vphase/2*self.vds)*self.fs
        
    def gate_f(self):
        ciss_0V = self.lsfet_params['Ciss_0V']
        qfet_gate = self.vgate*ciss_0V
        vbias = {'no':self.vgate,'yes':self.ckt_params['vin']}[self.ic_params['ldo']]
        return qfet_gate*self.vgate*(vbias/self.vgate)*self.fs
