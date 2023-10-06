
import pandas as pd
import pandas_flavor as pf

@pf.register_dataframe_method
def index_search(df,string):
    idx_str_search=df.index.str.lower().str.contains(string)
    return df.loc[idx_str_search]

@pf.register_dataframe_method
def index_exclude(df,string):
    idx_str_search=df.index.str.lower().str.contains(string)
    return df.loc[~idx_str_search]

@pf.register_dataframe_method
def column_search(df,column,string):
    idx_str_search=df[column].str.lower().str.contains(string)
    return df[idx_str_search]

@pf.register_dataframe_method
def column_exclude(df,column,string):
    idx_str_search=df[column].str.lower().str.contains(string)
    return df[~idx_str_search]

@pf.register_dataframe_method
def exclude_index_strings(df,string_list):
    """
    will exclude any factor provided in the string list
    """
    for str in string_list:
        df = df.index_exclude(str)
    return df

@pf.register_dataframe_method
def exclude_column_strings(df,column,string_list):
    """
    will exclude any factor provided in the string list
    """
    for str in string_list:
        df = df.column_exclude(column,str)
    return df

@pf.register_dataframe_method
def max_filter(df,column:str,maxval:float):
    """
    returns row with column value greater than maxval and sorts by that column
    """        
    return df[df[column]<maxval].sort_values(by=column)

@pf.register_dataframe_method
def min_filter(df,column:str,minval:float):
    """
    returns row with column value less than minval and sorts by that column
    """        
    return df[df[column]>minval].sort_values(by=column)