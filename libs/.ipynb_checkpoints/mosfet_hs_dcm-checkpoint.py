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

from unit_waveforms_turnon import Four_state


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
        self.cgd_0V = 420e-12 #self.fetparams['Ciss_0V']-self.c_gs() #may need to adjust this value if result of function is negative

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

class Fet_switching_on:   
    def __init__(self,args:tuple):  #idc and ipp are total current
        ic_params,hsfet_params,lsfet_params,vds,vgate,idc,ipp,m_hs,m_ls,rd = args
        self.icp = ic_params; self.hsfp = hsfet_params; self.lsfp = lsfet_params
        self.vds = vds; self.vdr = vgate
        self.idc = idc/m_hs; self.ipp = ipp/m_hs
        self.m_hs=m_hs; self.m_ls=m_ls; self.rd=rd
        
        self.vth = self.hsfp['Qgs']/self.hsfp['Ciss_Vds2'] 
        self.rghsdrvr = {5:self.icp['rg_hsdrvr_5V'],
                        10:self.icp['rg_hsdrvr_10V']}[self.vdr]
        self.rghs = self.hsfp['Rg']+self.m_hs*self.rghsdrvr    
        
        self.i_valley = max(0,self.idc-self.ipp/2); self.i_peak = max(self.ipp,self.idc+self.ipp/2)
        
        self.fet_cap = Fet_cap_vs_vds(self.hsfp,self.vds)
        #self.cgd = self.fet_cap.q_gd()/self.vds
        self.cds = self.fet_cap.c_ds(self.vds)
        self.cgs = self.fet_cap.c_gs()
        self.Cgd_0 = self.hsfp['Crss_Vds2']
        self.Cgd_1 = self.fet_cap.cgd_0V
        
        self.ld = self.hsfp['Ldrain']+(self.lsfp['Ldrain']+self.lsfp['Lsource'])/self.m_ls

        
    def period_td(self,**kwargs):  #function will be called twice - first for final state, then for f(t)
        tau = self.rghs*(self.Cgd_0+self.cgs)
        ton = -tau*ln(1-self.vth/self.vdr)
        if 't' in kwargs.keys():
            t=kwargs['t']
            return {'vgs_t':self.vdr*(1-e**-(t/tau)),
                    'id_t':0,
                    'vds_t':self.vds}
        else:
            return {'ton':ton}
        
    def period_t1(self,**kwargs):
        cgd = self.fet_cap.q_gd(self.vds)/self.vds
        gfs = self.hsfp['Gm']
        
        tau_gd = ((self.hsfp['Lsource']+self.ld)*self.rghs*cgd*gfs)**0.5
        tau_gs = gfs*self.hsfp['Lsource']+self.rghs*(cgd+self.cgs)
        exp_flag = (tau_gs > 2*tau_gd) #yes - exponential.  no - sinusoidal
        
        tau_a = 2*tau_gd**2/tau_gs
        w_a = ((1/tau_gd**2-(tau_gs/(2*tau_gd**2))**2)**2)**(1/4)   #sometimes would get imaginary
        tau_b = 2*tau_gd**2/(tau_gs-(tau_gs**2-4*(tau_gd)**2)**0.5)
        tau_c = 2*tau_gd**2/(tau_gs+(tau_gs**2-4*(tau_gd)**2)**0.5)

        def vgs_t_exp(t):
            a = tau_b*e**-(t/tau_b)-tau_c*e**-(t/tau_c)
            b = tau_b - tau_c
            return self.vdr-(self.vdr-self.vth)*a/b
        def vgs_t_sin(t):
            return self.vdr-(self.vdr-self.vth)*e**-(t/tau_a)*(cos(w_a*t)+sin(w_a*t)/w_a/tau_a)
        def id_t_exp(t):
            a = tau_b*e**-(t/tau_b)-tau_c*e**-(t/tau_c)
            b = tau_b - tau_c
            return gfs*(self.vdr-self.vth)*(1-a/b)
        def id_t_sin(t):
            return gfs*(self.vdr-self.vth)*(1-e**-(t/tau_a)*(cos(w_a*t)+sin(w_a*t)/w_a/tau_a))
        def vds_t_exp(t):
            a = e**-(t/tau_b)-e**-(t/tau_c)
            b = tau_b - tau_c
            return self.vds-gfs*(self.ld+self.hsfp['Lsource'])*(self.vdr-self.vth)*a/b
        def vds_t_sin(t):
            c = self.vds-gfs*(self.ld+self.hsfp['Lsource'])*(self.vdr-self.vth)*w_a*e**-(t/tau_a)
            d = (1+1/(w_a*tau_a)**2)*sin(w_a*t)
            return c*d
            
        def t_ir():
            def t_ir_max():
                est_ton = 0
                while id_t(est_ton) < self.i_valley:
                   est_ton += tau_gs/100
                return est_ton
            t_ir_ = t_ir_max()
            def t_ir_func(t):
                return id_t(t)-self.i_valley
            return opt.brentq(t_ir_func,0,t_ir_)  #root         
        def t_vf():
            def t_vf_max():
                est_ton = 0
                while (vds_t(est_ton) and (est_ton < ton_ir)):
                   est_ton += tau_gs/100
                return est_ton
            t_vf_ = t_vf_max()
            def t_vf_func(t):
                return vds_t(t)
            vds_flag = vds_t(t_vf_)>0
            if vds_flag:
                return t_vf_
            else:
                return opt.brentq(t_vf_func,0,t_vf_)

        #Below functions are f(t)
        vgs_t = {True:vgs_t_exp,False:vgs_t_sin}[exp_flag]
        id_t  = {True:id_t_exp, False:id_t_sin}[exp_flag]
        vds_t = {True:vds_t_exp,False:vds_t_sin}[exp_flag]
        
        if 't' in kwargs.keys():       
            t=kwargs['t']
            return {'vgs_t':vgs_t(t),
                    'id_t' :id_t(t),
                    'vds_t':vds_t(t)}                    
        else:   
            ton_ir = t_ir()
            ton_vf = t_vf()        
            i_slower_v = ton_vf < ton_ir
            ton_t1 = {True:ton_vf,False:ton_ir}[i_slower_v]
            id_on2_0 = {True:id_t(ton_t1),False:self.i_valley}[i_slower_v]
            vds_on2_0 = {True:0.1,False:vds_t(ton_t1)}[i_slower_v]
            vgs_on2_0 = vgs_t(ton_t1)

            return {'i_slower_v':i_slower_v,
                    'ton':ton_t1,
                    'id':id_on2_0,
                    'vds':vds_on2_0,
                    'vgs':vgs_on2_0}
        
    def period_t2(self,**kwargs): 
        gfs = self.hsfp['Gm']
        t0 = kwargs['init_conditions']
        cgd = self.fet_cap.q_gd(t0['vds'])/t0['vds'] 
        tau_gs = self.rghs*(cgd+self.cgs)
        ls = self.hsfp['Lsource']

        def vgs_t(t):
            return {True:(self.vdr-ls/(self.ld+ls)*self.vds-t0['vgs'])*(1-e**-(t/tau_gs))+t0['vgs'],
                    False:t0['vgs']}[t0['i_slower_v']]
        def id_t(t):
            return {True:t0['id']+self.vds/(self.ld+ls)*t,
                    False:self.i_valley}[t0['i_slower_v']]
        def vds_t(t):
            return {True:0,
                    False:t0['vds']-t*(gfs*(self.vdr-self.vth)-self.i_valley)/(1+gfs*self.rghs)/cgd}[t0['i_slower_v']]    
        def t_ir_max():
            est_ton = 0
            while id_t(est_ton) < self.i_valley:
               est_ton += tau_gs/100
            return est_ton
        def t_vf_max():
            est_ton = 0
            while (vds_t(est_ton) > 0):
               est_ton += tau_gs/100
            return est_ton
        def ton_t2():
            ir_max = t_ir_max()
            vf_max = t_vf_max()
            t2 = {True:  ir_max,
                  False: vf_max}[t0['i_slower_v']]
            def t_ir_func(t):
                return id_t(t)-self.i_valley
            def t_vf_func(t):
                return vds_t(t)
            return {True:  opt.brentq(t_ir_func,0,ir_max),
                    False: opt.brentq(t_vf_func,0,vf_max)}[t0['i_slower_v']]  #root

        if 't' in kwargs.keys():
            t=kwargs['t']
            return {'vgs_t':vgs_t(t),
                    'id_t':id_t(t),
                    'vds_t':vds_t(t)}
        else:
            t2 = ton_t2()
            id_on3_0  = id_t(t2)
            vds_on3_0 = vds_t(t2)
            vgs_on3_0 = vgs_t(t2)
            return {'i_slower_v':t0['i_slower_v'],
                    'ton':t2,
                    'id':id_on3_0,
                    'vds':vds_on3_0,
                    'vgs':vgs_on3_0}

    def period_t3(self,init_conditions:dict,t):
        t0 = init_conditions
        cgd = self.fet_cap.cgd_0V #fet_cap.q_gd(self.vds)/self.vds
        tau_gs = self.rghs*(cgd+self.cgs)
        def vgs_t(t):
            return (self.vdr-t0['vgs'])*(1-e**-(t/tau_gs))+t0['vgs']
        def id_t(t):
            return self.i_valley
        def vds_t(t):
            return 0
        return {'vgs_t':vgs_t(t),
                'id_t' :id_t(t),
                'vds_t':vds_t(t)}
        
    def whole_sequence(self,**kwargs): #t):
        #get final_states/init_conditions by calling functions WITHOUT providing 't' argument
        final_states = {'td':self.period_td(),
                        't1':self.period_t1()}
        final_states['t2']=self.period_t2(init_conditions=final_states['t1'])
        t_width  = {'td':final_states['td']['ton'],
                    't1':final_states['t1']['ton'],
                    't2':final_states['t2']['ton'],
                    't3':4e-8}
        self.t_widths = t_width
                
        if 't' in kwargs.keys():
            t=kwargs['t']
            pulses = Four_state(t_width).all_pulses(t) #use unit pulses for 4 periods - td,t1,t2,t3,t4
            p_td = self.period_td(t=t)
            p_t1 = self.period_t1(t=t-t_width['td'])#-p_td['ton'])
            p_t2 = self.period_t2(init_conditions=final_states['t1'],t=t-t_width['td']-t_width['t1'])
            p_t3 = self.period_t3(init_conditions=final_states['t2'],t=t-t_width['td']-t_width['t1']-t_width['t2'])
            p_t = {'td':p_td,'t1':p_t1,'t2':p_t2,'t3':p_t3}
             
            vgs_t = sum([p_t[tn]['vgs_t']*pulses[tn] for tn in p_t.keys()])
            id_t  = sum([p_t[tn]['id_t']*pulses[tn] for tn in p_t.keys()])
            vds_t = sum([p_t[tn]['vds_t']*pulses[tn] for tn in p_t.keys()])

            return {'vgs':vgs_t,
                    'id' :id_t,
                    'vds':vds_t}
        else:
            del t_width['t3']
            return t_width
    def e_on(self):
        tend = sum(self.whole_sequence().values())
        def p_on(t):
            vgs_id_vds_t = self.whole_sequence(t=t)
            return vgs_id_vds_t['vds']*vgs_id_vds_t['id']
        p_on_vectorized = np.vectorize(p_on)
        ti = np.linspace(0,tend,100)
        p_values = p_on_vectorized(ti)
        return np.trapz(p_values,ti)
    
    def plot_vgs_id_vds(self): #list of dataframes
        def whole_sequence_t(t):
            return self.whole_sequence(t=t)
        tfinal = 20e-9
        numdatapoints=1000  
        tstep = np.round(tfinal/numdatapoints,14)
        tarray=np.arange(0,tfinal,tstep, dtype=float)

        vg_id_vds=np.vectorize(whole_sequence_t)(tarray)
        vgs = [];id=[];vds=[];pds=[]
        for tslice in vg_id_vds.tolist():
            vgs.append(tslice['vgs'])
            id.append(tslice['id'])
            vds.append(tslice['vds'])
            pds.append(tslice['id']*tslice['vds']/10)
        vds[0]=self.vds

        f = plt.figure(figsize=(6,3))
        ax = f.add_subplot()
        
        ax.plot(tarray*1e9,vgs, 'r')
        ax.plot(tarray*1e9,id, 'b')
        ax.plot(tarray*1e9,vds, 'g')
        ax.plot(tarray*1e9,pds, 'm')

class Fet_switching_off:
    def __init__(self,args:tuple):      #idc and ipp are total current   
        ic_params,hsfet_params,lsfet_params,vds,vgate,idc,ipp,m_hs,m_ls,rd = args
        self.icp = ic_params; self.hsfp = hsfet_params; self.lsfp = lsfet_params        
        self.vds = vds; self.vdr = vgate        
        self.idc = idc/m_hs; self.ipp = ipp/m_hs
        self.m_hs=m_hs;self.m_ls=m_ls
        self.rd = rd
        
        self.rghsdrvr = {5:self.icp['rg_hsdrvr_5V'],
                        10:self.icp['rg_hsdrvr_10V']}[self.vdr]
        self.rghs = self.hsfp['Rg']+self.m_hs*self.rghsdrvr        
        self.vth = self.hsfp['Qgs']/self.hsfp['Ciss_Vds2']
        self.i_valley = self.idc-self.ipp/2; self.i_peak = self.idc+self.ipp/2
        
        self.fet_cap = Fet_cap_vs_vds(self.hsfp,self.vds)
        #self.cgd = self.fet_cap.q_gd()/self.vds
        self.cds = self.fet_cap.c_ds(self.vds)
        self.cgs = self.fet_cap.c_gs()
        self.Cgd_0 = self.hsfp['Crss_Vds2']
        self.Cgd_1 = self.fet_cap.cgd_0V
        
        self.ld = self.hsfp['Ldrain']+(self.lsfp['Ldrain']+self.lsfp['Lsource'])/self.m_ls

    def period_td(self,**kwargs):
        tau_gs = self.rghs*(self.Cgd_1+self.cgs)
        gfs = self.hsfp['Gm']
        
        def vgs_t(t):
            return self.vdr*e**-(t/tau_gs)
        id_t = self.i_peak
        vds_t = 0

        def t_td_max():            
            est_t = 0
            while vgs_t(est_t) > self.vth + self.i_peak/gfs:
               est_t += tau_gs/100
            return est_t
        def t_td_func(t):
            return vgs_t(t)-self.vth-self.i_peak/gfs
        t_td = opt.brentq(t_td_func,0,t_td_max())

        if 't' in kwargs.keys():
            t=kwargs['t']
            return {'vgs_t':vgs_t(t),
                    'id_t': id_t,
                    'vds_t':vds_t}
        else:
            return {'ton':t_td,
                    'vgs':vgs_t(t_td),
                    'id': id_t,
                    'vds':vds_t}

    def period_t1(self,**kwargs):
        gfs = self.hsfp['Gm']
        t0 = kwargs['init_conditions']
        cgd = self.fet_cap.q_gd(self.vds)/self.vds
        tau_gs = self.rghs*(cgd+self.cgs)
        
        vgs_t = t0['vgs']
        id_t = self.i_peak
        def vds_t(t):
            a=gfs*self.vth+self.i_peak
            b=(1+gfs*self.rghs)*cgd
            return a/b*t

        def t_t1_max():            
            est_t = 0
            while vds_t(est_t) <self.vds:
               est_t += tau_gs/100
            return est_t
        def t_t1_func(t):
            return vds_t(t)-self.vds
        t_t1 = opt.brentq(t_t1_func,0,t_t1_max())

        if 't' in kwargs.keys():
            t=kwargs['t']
            return {'vgs_t':vgs_t,
                    'id_t': id_t,
                    'vds_t':vds_t(t)}
        else:
            return {'ton':t_t1,
                    'vgs':vgs_t,
                    'id': id_t,
                    'vds':vds_t(t_t1)}

    def period_t2(self,**kwargs):
        gfs = self.hsfp['Gm']
        t0 = kwargs['init_conditions']
        cgd = self.fet_cap.q_gd(self.vds)/self.vds
        tau_gd = ((self.hsfp['Lsource']+self.ld)*self.rghs*self.Cgd_0*gfs)**0.5
        tau_gs = gfs*self.hsfp['Lsource']+self.rghs*(self.Cgd_0+self.cgs)
        exp_flag = (tau_gs > 2*tau_gd) #yes - exponential.  no - sinusoidal
        
        tau_a = 2*tau_gd**2/tau_gs
        w_a = ((1/tau_gd**2-(tau_gs/(2*tau_gd**2))**2)**2)**(1/4)
        tau_b = 2*tau_gd**2/(tau_gs-(tau_gs**2-4*(tau_gd)**2)**0.5)
        tau_c = 2*tau_gd**2/(tau_gs+(tau_gs**2-4*(tau_gd)**2)**0.5)

        def vgs_t_exp(t):
            a = tau_b*e**-(t/tau_b)-tau_c*e**-(t/tau_c)
            b = tau_b - tau_c
            return (self.i_peak/gfs+self.vth)*a/b
        def vgs_t_sin(t):
            a = (self.i_peak/gfs+self.vth)*e**-(t/tau_a)
            return a*(cos(w_a*t)+sin(w_a*t)/w_a/tau_a)
            
        def id_t_exp(t):
            a = tau_b*e**-(t/tau_b)-tau_c*e**-(t/tau_c)
            b = tau_b - tau_c
            return (gfs*self.vth+self.i_peak)*(a/b)-gfs*self.vth
        def id_t_sin(t):
            a=(gfs*self.vth+self.i_peak)*e**-(t/tau_a)
            return a*(cos(w_a*t)+sin(w_a*t)/w_a/tau_a)
            
        def vds_t_exp(t):
            a = e**-(t/tau_b)-e**-(t/tau_c)
            b = tau_b - tau_c
            return self.vds+(self.ld+self.hsfp['Lsource'])*(gfs*self.vth+self.i_peak)*a/b
        def vds_t_sin(t):
            c = w_a*e**-(t/tau_a)
            d = (1+1/(w_a*tau_a)**2)*sin(w_a*t)
            return self.vds+(self.ld+self.hsfp['Lsource'])*(gfs*self.vth+self.i_peak)*c*d
        #Below functions are f(t)
        vgs_t = {True:vgs_t_exp,False:vgs_t_sin}[exp_flag]
        id_t  = {True:id_t_exp, False:id_t_sin}[exp_flag]
        vds_t = {True:vds_t_exp,False:vds_t_sin}[exp_flag]

        def t_t2_max():
            est_t = 0
            while id_t(est_t) >0:
               est_t += tau_gs/100
            return est_t
        def t_t2_func(t):
            return id_t(t)
        t_t2 = opt.brentq(t_t2_func,0,t_t2_max())

        self.vds_ringpk = vds_t(t_t2)
             
        if 't' in kwargs.keys():
            t=kwargs['t']
            return {'vgs_t':vgs_t(t),
                    'id_t': id_t(t),
                    'vds_t':vds_t(t)}
        else:
            return {'ton':t_t2,
                    'vgs':vgs_t(t_t2),
                    'vds_ringpk':self.vds_ringpk,
                    'vds':vds_t(t_t2)}

    def period_t3(self,init_conditions,t):
        t0=init_conditions
        cgd = self.hsfp['Crss_Vds2']
        cds = self.cds
        tau_gs = self.rghs*(cgd+self.cgs)
        tau_D = 2*self.ld/self.rd
        w_D=(1/(self.ld*(cgd+cds))-self.rd/2/self.ld)**0.5
             
        vgs_t=t0['vgs']*e**-(t/tau_gs)
        def vds_t_func(t):
            return self.vds+(t0['vds']-self.vds)*e**-(t/tau_D)*cos(w_D*t)
        vds_t = vds_t_func(t)
        id_t = (cgd+cds)*nd.Derivative(vds_t_func,step=1e-10)(t)
        return {'cgd':cgd,
                'cgs':self.cgs,
                'vgs_t':vgs_t,
                'id_t':id_t,
                'vds_t':vds_t}

    def whole_sequence(self,**kwargs): #t):
        #get final_states/init_conditions by calling functions WITHOUT providing 't' argument
        final_states = {'td':self.period_td()}
        final_states['t1'] = self.period_t1(init_conditions=final_states['td'])
        final_states['t2'] = self.period_t2(init_conditions=final_states['t1'])
        t_width  = {'td':final_states['td']['ton'],
                    't1':final_states['t1']['ton'],
                    't2':final_states['t2']['ton'],
                    't3':4e-8}
        self.t_widths = t_width

        if 't' in kwargs.keys():
            t=kwargs['t']
            pulses = Four_state(t_width).all_pulses(t) #use unit pulses for 4 periods - td,t1,t2,t3,t4
    
            p_td = self.period_td(t=t)
            p_t1 = self.period_t1(init_conditions=final_states['td'],t=t-t_width['td'])#-p_td['ton'])
            p_t2 = self.period_t2(init_conditions=final_states['t1'],t=t-t_width['td']-t_width['t1'])
            p_t3 = self.period_t3(init_conditions=final_states['t2'],t=t-t_width['td']-t_width['t1']-t_width['t2'])
            p_t = {'td':p_td,'t1':p_t1,'t2':p_t2,'t3':p_t3}
            self.periods_t = p_t
             
            vgs_t = sum([p_t[tn]['vgs_t']*pulses[tn] for tn in p_t.keys()])
            id_t  = sum([p_t[tn]['id_t']*pulses[tn] for tn in p_t.keys()])
            vds_t = sum([p_t[tn]['vds_t']*pulses[tn] for tn in p_t.keys()])
    
            return {'vgs':vgs_t,
                    'id' :id_t,
                    'vds':vds_t}
        else:
            del t_width['t3']
            return t_width

    def e_off(self):
        tend = sum(self.whole_sequence().values())
        def p_off(t):
            vgs_id_vds_t = self.whole_sequence(t=t)
            return vgs_id_vds_t['vds']*vgs_id_vds_t['id']
        p_off_vectorized = np.vectorize(p_off)
        ti = np.linspace(0,tend,100)
        p_values = p_off_vectorized(ti)
        return np.trapz(p_values,ti)

    def plot_vgs_id_vds(self): #list of dataframes
        def whole_sequence_t(t):
            return self.whole_sequence(t=t)
        tfinal = 40e-9
        numdatapoints=1000  
        tstep = np.round(tfinal/numdatapoints,14)
        tarray=np.arange(0,tfinal,tstep, dtype=float)

        vg_id_vds=np.vectorize(whole_sequence_t)(tarray)
        vgs = [];id=[];vds=[];pds=[]
        for tslice in vg_id_vds.tolist():
            vgs.append(tslice['vgs'])
            id.append(tslice['id'])
            vds.append(tslice['vds'])
            pds.append(tslice['id']*tslice['vds']/10)
        vgs[0]=vgs[1]
        id[0]=id[1]

        f = plt.figure(figsize=(6,3))
        ax = f.add_subplot()
        
        ax.plot(tarray*1e9,vgs, 'r')
        ax.plot(tarray*1e9,id, 'b')
        ax.plot(tarray*1e9,vds, 'g')
        ax.plot(tarray*1e9,pds, 'm')

class Losses:
    def __init__(self,ckt_params,fs_dcm,*args): 
        self.args=args
        self.ic_params,self.hsfet_params,lsfet_params,self.vds,self.vgate,self.idc,self.ipp,self.m_hs,m_ls,rd = self.args        
        self.ckt_params = ckt_params 
        self.state_count = self.ckt_params['state count']
        self.ts = {4:(2*self.ckt_params['t_state13']+2*self.ckt_params['t_state24']),
                   2:self.ckt_params['t_state13']+self.ckt_params['t_state24']}[self.state_count]
        self.fs=fs_dcm*{2:1,4:0.5}[self.state_count]   #1/self.ts        
        self.fet_switch_on_obj = Fet_switching_on(self.args)
        self.fet_switch_off_obj = Fet_switching_off(self.args)        
        self.summary = {'sw_on':self.fs*self.fet_switch_on_obj.e_on(),
                        'sw_off':self.fs*self.fet_switch_off_obj.e_off(),
                        'cond': self.cond_f(),
                        'ring': self.ring_f(),
                        'gate': self.gate_f()}
    
    def cond_f(self):
        tcoeff = 3500e-6
        tmult = tcoeff*(self.ckt_params['Tamb']-25)
        rdson = {5:self.hsfet_params['Rdson_4.5V'],10:self.hsfet_params['Rdson_10V']}[self.vgate]*(1+tmult)
        t_Qhs = self.ckt_params['t_Qhs']
        i_fetrms = (((self.idc**2+self.ipp**2/12)*t_Qhs/self.ts)**0.5)/self.m_hs
        return i_fetrms**2*rdson 

    def ring_f(self):
        fet_ds_ring_pk = self.fet_switch_off_obj.vds_ringpk #this must be called after e_off
        qoss_vphase = Fet_cap_vs_vds(self.hsfet_params,self.vds).q_oss(self.vds)
        qoss_vds_ring_pk = Fet_cap_vs_vds(self.hsfet_params,fet_ds_ring_pk).q_oss(fet_ds_ring_pk)        
        return (qoss_vds_ring_pk/2*(fet_ds_ring_pk-2*self.vds)+qoss_vphase/2*self.vds)*self.fs
        
    def gate_f(self):
        ciss_0V = self.hsfet_params['Ciss_0V']
        qfet_gate = self.vgate*ciss_0V
        vbias = {'no':self.vgate,'yes':self.ckt_params['vin']}[self.ic_params['ldo']]
        return qfet_gate*self.vgate*(1/2+vbias/self.vgate/2)*self.fs
