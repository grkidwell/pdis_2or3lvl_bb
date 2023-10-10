
from circuit_4state import circuit_params as circuit_params_4state
from circuit_2state import circuit_params as circuit_params_2state

from controller import get_ic_params
from cyntec import Inductor_pdis, create_ind_family_df, ind_pdis_obj, l_set, cyntec_filename
import mosfet_hs_dcm as mosfet_hs
import mosfet_ls_dcm as mosfet_ls # import get_fet_params, Fet_cap_vs_vds, Fet_switching_on, Fet_switching_off, Losses
import mosfet_q4_dcm as mosfet_q4


from unit_waveforms_turnon import Four_state

from capacitors import Caplosses

from board import p_cu
from current_shunt import p_shunt


class Buckconverter_losses:
    def __init__(self,inp_params:dict):
        self.ip = inp_params
        self.cont_p = get_ic_params(self.ip['controller'])
        self.lvl_config = self.ip['lvl_config']        
        self.ckt_params = {'2 level':circuit_params_2state(self.ip['vin'],self.ip['vout'],self.ip['iout'],self.ip['fs'],
                                                           self.ip['tambient'],self.ip['lout']['config']),
                           '3 level':circuit_params_4state(self.ip['vin'],self.ip['vout'],self.ip['iout'],self.ip['fs'],
                                                           self.ip['tambient'],self.ip['lout']['config'])}\
                           [self.lvl_config]
        self.vds = self.ckt_params['vphase']
        self.vgate = self.ip['vgate']
        self.hs_p = mosfet_hs.get_fet_params(self.ip['hsfet_partnum'])
        self.ls_p = mosfet_ls.get_fet_params(self.ip['lsfet_partnum'])
        self.q4_p = mosfet_q4.get_fet_params(self.ip['q4_partnum'])
        
        self.lout_obj = ind_pdis_obj(self.ckt_params,self.ip['lout']['value(uH)'],create_ind_family_df(self.ip['lout']['family']))
        self.fs_dcm = self.lout_obj.fs_dcm   #inductor ripple frequency
        self.idc = self.lout_obj.idc*{'single':1,'series':1,'parallel':2}[self.ip['lout']['config']]
        self.ipp = self.lout_obj.ipp*{'single':1,'series':1,'parallel':2}[self.ip['lout']['config']]

        self.p_lout = self.lout_obj.summary
        self.hs_Losses_obj = mosfet_hs.Losses(self.ckt_params,self.fs_dcm,self.cont_p,self.hs_p,self.ls_p,self.vds,self.vgate,self.idc,self.ipp,
                                     self.ip['m_hs'],self.ip['m_ls'],self.ip['rd'])
        self.p_hs = self.hs_Losses_obj.summary 
        self.ls_losses_obj = mosfet_ls.Losses(self.ckt_params,self.fs_dcm,self.cont_p,self.hs_p,self.ls_p,self.vds,self.vgate,self.idc,self.ipp,
                                     self.ip['m_hs'],self.ip['m_ls'],self.ip['rd'])
        self.p_ls = self.ls_losses_obj.summary
        self.q4_Losses_obj = mosfet_q4.Losses(self.ip,self.q4_p,self.lout_obj.irms_dcm)
        self.p_q4 = self.q4_Losses_obj.summary

        
        #caploss and cu calcs still need to be modified to account for period stretching
        self.p_caps = Caplosses(self.ckt_params,self.idc,self.ipp,self.ip['caps']).summary        
        self.p_cu = p_cu(self.ckt_params['Iinp'],self.ckt_params['Idc'],self.ip['tambient'])
        self.p_ic = self.cont_p['pic']

        self.p_summary ={'hs turn-on':self.p_hs['sw_on'],
                         'hs turn-off':self.p_hs['sw_off'],
                         'hs rdson':self.p_hs['cond'],
                         'hs ringing':self.p_hs['ring'],
                         'hs gate':self.p_hs['gate'],
                         'ls rdson':self.p_ls['cond'],
                         'ls bd':self.p_ls['bd_on']+self.p_ls['bd_off'],
                         'ls ring_qrr':self.p_ls['ring'],
                         'ls gate':self.p_ls['gate'],
                         'q4_rdson':self.p_q4['cond'],
                         'lout rdc+rac':self.p_lout['dcr'],
                         'lout core':self.p_lout['core'],
                         'flying cap':self.p_caps['flying'],
                         'input cap':self.p_caps['vin'],
                         'board cu':self.p_cu,
                         'ic':self.p_ic}
        self.p_summary = {key:round(value,3) for key,value in self.p_summary.items()}

        self.p_totals ={'hs fet':sum(list(self.p_hs.values()))-self.p_summary['hs gate'],
                        'ls fet':sum(list(self.p_ls.values()))-self.p_summary['ls gate'],
                        'q4 fet':sum(list(self.p_q4.values())),
                        'lout':sum(list(self.p_lout.values())),
                        'caps':sum(list(self.p_caps.values())),
                        'ic_with_gate' : self.p_summary['ic']+(self.p_summary['hs gate']*self.ip['m_hs']+self.p_summary['ls gate']*self.ip['m_ls'])}

        
        self.ptotal = (self.p_totals['hs fet']*self.ip['m_hs']+self.p_totals['ls fet']*self.ip['m_ls'])*{'2 level':1,'3 level':2}[self.lvl_config]+ \
                  self.p_totals['q4 fet']+\
                  self.p_totals['lout']*(2)**(self.ip['lout']['config']!='single')+\
                  self.p_totals['caps']+self.p_summary['board cu']+self.p_totals['ic_with_gate']

        #optional components like shunts:
        self.iin = self.ip['vout']*self.ip['iout']/self.ip['vin']
        self.input_shunt_calc() 

        
        self.efficiency = round(self.ip['vout']*self.ip['iout']/(self.ip['vout']*self.ip['iout'] + self.ptotal),3)
        self.p_totals['total']=self.ptotal
        self.p_totals['efficiency']=self.efficiency
        self.p_totals = {key:round(value,3) for key,value in self.p_totals.items()}
                
    def input_shunt_calc(self):
        if 'r_shunt_input' in self.ip.keys():
            rshunt = self.ip['r_shunt_input']
            pshunt = p_shunt(rshunt,self.iin)
            self.p_summary['inp_shunt']=pshunt
            self.ptotal = self.ptotal+pshunt
        
    
def b_obj_var(inp_params:dict,variable:dict):
    p = inp_params
    p[variable['param']] = variable['value']
    return Buckconverter_losses(p) 

def df_var(inp_params,variable:dict):
    b_obj = b_obj_var(inp_params,variable)
    return pd.DataFrame.from_dict(b_obj.p_summary| b_obj.p_totals,orient='index',columns=[variable['value']]).T

def df_var_list(inp_params:dict,varlist:dict):
    param = varlist['param']
    value_list = varlist['values']
    return [df_var(inp_params,{'param':param,'value':value}) for value in value_list]



