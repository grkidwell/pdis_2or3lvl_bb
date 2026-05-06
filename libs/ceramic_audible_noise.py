
#these functions operate in dataframe created by capacitor_selection_tables.py library
import math

def ceq(n,series,cap):
    if series and n%2:
        print("series config quantities must be multiples of 2")
    cdiv = 2**series
    cgroups = n/cdiv
    return (cap/cdiv)*cgroups

# def dv_factor(n,series,cap):
#     c_eq = ceq(n,series,cap)
#     return 1/(c_eq*2**series)

# def dv_factor(series):
#     return 1/(2**series)

def dx(n,series,package):
    if package in ['B','D','X','V']:
        w = {'B':11,'D':17,'X':17,'V':17}[package]
    else:
        w = int(str(package)[2:]) 
    cdiv = 2**series
    cgroups = n/cdiv
    return w*cgroups

def dy(series,package):
    if package in ['B','D','X','V']:
        length = {'B':13,'D':29,'X':29,'V':29}[package]
    else:
        length = int(str(package)[0:2]) 
    cdiv = 2**series
    return length*cdiv

aud_noise_factor_db = {'0402':-6,'0603':0,'0805':6.5,'1206':12.5,'1210':20}   #per cap

def db_to_energy_ratio(db):
    return 10**(db/10)

def energy_ratio_to_db(e):
    return math.log10(e)*10

def vpp_ratio_to_db(vppratio):    #vpp
    return math.log10(vppratio)*20

def sound_factor_db(n,series,db_0603,vppratio,package,partnum):
    if 'K' in partnum[0]:
        low_noise_adj = -20
    elif 'Z' in partnum[0]:
        low_noise_adj = -15
    else:
        low_noise_adj = 0
    if package not in ['B','D','X','V']:
        vpp_ratio = vppratio/2**series
        db = db_0603+vpp_ratio_to_db(vpp_ratio)+aud_noise_factor_db[package]+low_noise_adj
        e = db_to_energy_ratio(db)
        etotal = n*e
        er2db = round(energy_ratio_to_db(etotal),1)
    else:
        er2db = 0
    return er2db

def add_columns(df):
    df['ceq'] = df.apply(lambda row: ceq(row['qnty'],row['series'],row['cap_uF(@vbias)']),axis=1)
    #df['dv_factor']= df.apply(lambda row: dv_factor(row['series']),axis=1)
    df['x'] = df.apply(lambda row: dx(row['qnty'],row['series'],row['package']),axis=1)
    df['y'] = df.apply(lambda row: dy(row['series'],row['package']),axis=1)

def add_sound_factor(df,db_0603,vppratio):
    df['sound factor(dB)'] = df.apply(lambda row: sound_factor_db(row['qnty'],row['series'],db_0603,vppratio,row['package'],row['partnum']),axis=1)