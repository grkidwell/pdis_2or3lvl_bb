import pandas as pd
import numpy as np
import pathlib

def create_datafilelist():
    path = r'data/murata_capacitors'
    datafiles = [file for file in list(pathlib.Path(path).iterdir()) if 'csv' in file.name]
    return datafiles

def getpartnum(posixpath):
    return posixpath.name.split('_')[0] 
    
def create_df_of_folder_contents(datafilelist):
    folder_content_list = [getpartnum(datafile) for datafile in datafilelist]
    folder_content_dict = {partnum:murata_pn_parser(partnum) for partnum in folder_content_list}
    df2 = pd.DataFrame.from_dict(folder_content_dict,orient='index').reset_index().rename(columns={'index':'partnum'})
    df2['manufacturer']='murata'
    df2['type']='mlcc'
    return df2[['manufacturer','partnum','type','package','temp','volt rating','cap_uF(@0V)']]
    

def murata_pn_parser(pn):
    pn_dict1   = {#'series':pn[0:3],
                  'package':pn[3:5],
                  'temp':pn[6:8],
                  'volt rating':pn[8:10],
                  'cap_uF(@0V)':pn[10:13]}
    decoder_dict={
    #'series'  : {pn_dict1['series']:pn_dict1['series']},   
    #'type':
    'package' : {'15':'0402',
                    '18':'0603',
                    '21':'0805',
                    '31':'1206',
                    '32':'1210'},
    'temp'    : {'R6':'X5R',
                 'R7':'X7R',
                 'C7':'X7S',
                 'C8':'X6S',
                 'D7':'X7T',
                 'D8':'X6T',
                 'Z7':'X7R'},
    'volt rating' : {'OJ':6.3,
                 '1A':10,
                 '1C':16,
                 '1E':25,
                 'YA':35,
                 '1H':50,
                 '2A':100},
    'cap_uF(@0V)'  : {pn_dict1['cap_uF(@0V)']:pn_dict1['cap_uF(@0V)']}}
    pn_dict2={param:decoder_dict[param][value] for param,value in pn_dict1.items()}
    return pn_dict2
    
def cap_in_uF(strval):
    val = float(strval[:2])
    exp = float(strval[2])
    return val*10**(exp-6) 
    

class Ceramic_cap:
    def __init__(self,datafile,vbias):
        self.partnum = getpartnum(datafile)
        self.vbias = vbias
        df = pd.read_csv(datafile,skiprows=4)
        self.df_cap_vs_bias = self.create_df_cap_vs_bias(df)
        self.df_param = self.create_param_df(df)
        #self.pn_parser = murata_pn_parser(self.partnum)

    def create_df_cap_vs_bias(self,df):
        df2=df.copy()
        df2.columns = ['DCbias(V)','Capacitance(uF)','duh']
        df2 = df2.iloc[3:].drop(columns=['duh'])
        df2['DCbias(V)']=df2['DCbias(V)'].astype(float)
        df2['Capacitance(uF)']=df2['Capacitance(uF)'].astype(float)*1e6
        return df2
        
    def get_derated_cap_value(self,df,voltage): 
        vdelta = np.abs(df['DCbias(V)'] - voltage)
        nearest_index = vdelta.idxmin()
        return round(float(df['Capacitance(uF)'].iloc[nearest_index]),3)

    def create_param_df(self,df):
        df_params = df.copy().iloc[0:2,0:2]
        df_params = df_params.T.reset_index(drop=True)
        df_params.columns = df_params.iloc[0]
        df_params = df_params.drop(index=0).reset_index(drop=True)
        df_params['partnum']=self.partnum
        pn_parser = murata_pn_parser(self.partnum)
        df_pnparse = pd.DataFrame.from_dict(pn_parser,orient='index').T
        df_params = df_params.join(df_pnparse)
        df_params['cap_uF(@vbias)'] = self.get_derated_cap_value(self.df_cap_vs_bias,self.vbias)
        df_params['cap_uF(@0V)']=cap_in_uF(df_params['cap_uF(@0V)'].tolist()[0])
        return df_params[['partnum','package','temp','volt rating','cap_uF(@0V)','cap_uF(@vbias)','ESR_mOhm','ESL_nH']]

    def plot_cap_derating(self):
        df = self.df_cap_vs_bias
        ax = df.plot(x='DCbias(V)',y='Capacitance(uF)',title='Cap Derating',legend=False)
        ax.set_ylabel('Capacitance(uF)')
        ax.grid(True)
        


