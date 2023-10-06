import pandas as pd

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


class Losses:
    def __init__(self,input_params,q3_params,lout_irms_dcm):
        self.inp_params = input_params.copy()
        self.q3_params = q3_params.copy()
        self.vgate = self.inp_params['vgate']
        lout_config = self.inp_params['lout']['config']
        self.i_fetrms = lout_irms_dcm*{'single':1,'series':1,'parallel':2}[lout_config]
        self.summary = {'cond':self.cond_f()}

    def cond_f(self):
        tcoeff = 3500e-6
        tmult = tcoeff*(self.inp_params['tambient']-25)
        rdson = {5:self.q3_params['Rdson_4.5V'],10:self.q3_params['Rdson_10V']}[self.vgate]*(1+tmult)
        return self.i_fetrms**2*rdson 
        