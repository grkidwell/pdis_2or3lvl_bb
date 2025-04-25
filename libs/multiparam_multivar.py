import libs.append_path
from add_python_libraries import *

from bucklosses_dcm import Buckconverter_losses

class Multiparam_multivar_df:
    def __init__(self,input_params:dict,inp_config:dict,varlist:dict):
        self.input_params = input_params
        self.inp_config = inp_config
        self.varlist = varlist
        self.df = self.create_df()
        
    def b_obj_var(self,inp_params:dict,variable:dict):
        p = inp_params
        p[variable['param']] = variable['value']
        return Buckconverter_losses(p) 
    
    def df_var(self,inp_params,variable:dict):
        b_obj = self.b_obj_var(inp_params,variable)
        return pd.DataFrame.from_dict(b_obj.p_summary| b_obj.p_totals,orient='index',columns=[variable['value']]).T
    
    def df_var_list(self,inp_params:dict,varlist:dict):
        param = varlist['param']
        value_list = varlist['values']
        return [self.df_var(inp_params,{'param':param,'value':value}) for value in value_list]
    
    def mod_inp_params(self,inp_params,subconfig):
        new_config = inp_params.copy()
        for key,value in subconfig.items():
            new_config[key]=value
        return new_config

    def create_df(self):
        new_params={config_idx:self.mod_inp_params(self.input_params,configx) for config_idx,configx in self.inp_config.items()}
        df_pdis_results = {config_idx:pd.concat(self.df_var_list(new_params_x,self.varlist)) for config_idx,new_params_x in new_params.items()}
        for config_idx,df in df_pdis_results.items():
            df['config']=f'config{config_idx}'
        df_all = pd.concat(df_pdis_results.values())
        return df_all

    def parametric_plot_eff(self):
        self.df.drop([0]).groupby('config')['efficiency'].plot(legend=True)

    def parametric_plot_pdis(self):
        self.df.drop([0]).groupby('config')['total'].plot(legend=True)

    def parametric_plot_Fs(self):
        self.df.drop([0]).groupby('config')['Fs'].plot(legend=True)