import pandas as pd

ic_filename = r'data/ic_params.xlsx'

def get_ic_params(partnumber:str):
    df_ic = pd.read_excel(ic_filename)
    df = df_ic
    paramdict = dict(df[['parameter',partnumber]].values)
    unitdict = dict(df[['parameter','units']].values)
    scalefactors = dict(df[['parameter','multiplier']].values)
    params = {param:value*scalefactors[param] for param,value in paramdict.items() if unitdict[param] != 'na'}
    params['package']=paramdict['package']
    params['ldo']=paramdict['ldo']
    return params