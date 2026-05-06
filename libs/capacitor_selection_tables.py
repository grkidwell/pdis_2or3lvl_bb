
import pandas as pd
from ipydatagrid import DataGrid
from json import load

from ceramic_cap import create_datafilelist, create_df_of_folder_contents, Ceramic_cap

class Cap_selector:
    def __init__(self,vbias):
        self.vbias = vbias
        self.datafiles = create_datafilelist()

    def create_mlcc_selection_table(self):   
        df = create_df_of_folder_contents(self.datafiles)
        df['qnty']=0; df['series']=0
        df=df[['manufacturer','partnum','qnty','series','type','package','temp','volt rating','cap_uF(@0V)']]
        return DataGrid(df, editable=True, auto_fit_columns=True,layout={"height":"500px", "width":"1000px"})

    def create_mlcc_selection_table_ab(self):   
        df = create_df_of_folder_contents(self.datafiles)
        df['qnty a']=0; df['series a']=0
        df['qnty b']=0; df['series b']=0
        df=df[['manufacturer','partnum','qnty a','series a','qnty b','series b','type','package','temp','volt rating','cap_uF(@0V)']]
        return DataGrid(df, editable=True, auto_fit_columns=True,layout={"height":"500px", "width":"1200px"})

    def create_sp_selection_table(self):
        filename = r'data/SP_caps/SP_captable.xlsx'
        df_spcaps = pd.read_excel(filename)
        df_spcaps['qnty']=0#; df_spcaps['qnty b']=0
        df_spcaps=df_spcaps[['manufacturer','partnum','qnty','type','package','temp','volt rating','cap_uF(@vbias)','ESR_mOhm','ESL_nH','Irms_A']]
        return DataGrid(df_spcaps, editable=True, auto_fit_columns=True,layout={"height":"350px", "width":"1200px"})

    def create_sp_selection_table_ab(self):
        filename = r'data/SP_caps/SP_captable.xlsx'
        df_spcaps = pd.read_excel(filename)
        df_spcaps['qnty a']=0; df_spcaps['qnty b']=0;df_spcaps['series a']=0; df_spcaps['series b']=0
        df_spcaps=df_spcaps[['manufacturer','partnum','qnty a','qnty b', 'type','package','temp','volt rating','cap_uF(@vbias)','ESR_mOhm','ESL_nH','Irms_A','series a','series b']]
        return DataGrid(df_spcaps, editable=True, auto_fit_columns=True,layout={"height":"350px", "width":"1200px"})

    def create_allcaps_selected_df(self,mlcc_datagrid,sp_datagrid,a_or_b:str):
        def df_a_or_b(df,a_or_b:str):
            q_old = f'qnty {a_or_b}'
            s_old = f'series {a_or_b}'
            if f'series {a_or_b}' in df.columns.tolist():
                df_new = df[['key','partnum','type',q_old,s_old]][df[q_old]>0].rename(columns={q_old:'qnty',s_old:'series'})
            else:
                df_new = df[['key','partnum','type',q_old]][df[q_old]>0].rename(columns={q_old:'qnty'})
            return df_new
        def create_mlcc_selected_df(mlcc_datagrid):
            df2=df_a_or_b(mlcc_datagrid.data.reset_index(),a_or_b)
            # if ab:
            #     df_caplist = df2[(df2['qnty a']>0) | (df2['qnty b']>0)].reset_index()[['key','partnum','type','qnty a','series a','qnty b','series b']]
            # else:
            df_caplist = df2[df2.qnty>0][['key','partnum','type','qnty','series']]
            self.cap_obj_list = [Ceramic_cap(self.datafiles[df_cap.key],self.vbias/2**df_cap.series) for _,df_cap in df_caplist.iterrows()]
            df_duh = pd.concat([capobj.df_param for capobj in self.cap_obj_list]).reset_index(drop=True)
            df_mlcc_selected = pd.merge(df_caplist,df_duh,on='partnum')
            df_mlcc_selected[['volt rating','ESR_mOhm','ESL_nH']]=df_mlcc_selected[['volt rating','ESR_mOhm','ESL_nH']].astype(float)
            return df_mlcc_selected
        def create_sp_selected_df(sp_datagrid):
            df2=sp_datagrid.data.reset_index()
            df3=df_a_or_b(df2,a_or_b)
            df_sp_selected = pd.merge(df2.drop(['manufacturer','qnty a','qnty b','type','key'],axis=1),df3[df3.qnty>0],on='partnum') 
            return df_sp_selected.drop(columns = ['series a','series b'])
        df_mlcc_selected=create_mlcc_selected_df(mlcc_datagrid)
        df_sp_selected=create_sp_selected_df(sp_datagrid)
        #with pd.option_context('future.no_silent_downcating', True):
        df = pd.concat([df_mlcc_selected,df_sp_selected])
        df_mask = df.mask(df.isna(),'') #fillna(0)
        return df_mask.reset_index(drop=True)


    def save_to_file(self,dg_obj,filename):
        dg_obj.data.to_csv(filename)

    def load_from_file(self,dg_obj,filename):
        df = pd.read_csv(filename,index_col=False)
        df.drop('key',axis=1, inplace=True)
        dg_obj.data = df
        
