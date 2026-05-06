import libs.append_path
from add_python_libraries import *
#import pd_filter_fcns
import pandas as pd

import os, datetime
from openpyxl import load_workbook
from openpyxl.drawing.image import Image
from io import BytesIO

from capacitorZ_multibank import Zcapbanks, Zcapbanks_df

from current_buck_fft import Current_buck_FFT

from current_capbanks_fft import Current_capbanks_FFT
from current_capbanks_waveforms import Current_capbanks_waveforms

from voltage_capbanks_waveform import Vripple_capbanks_waveforms

from ceramic_audible_noise import db_to_energy_ratio,energy_ratio_to_db,add_columns,add_sound_factor

class InputCapbanks:
    def __init__(self,iobj,df_allcapsselected):
        self.iobj = iobj
        self.currents_buckfftobj = Current_buck_FFT(self.iobj)
        self.df_allcaps_selected = df_allcapsselected.copy()
        add_columns(self.df_allcaps_selected)
        self.zcapbanks_obj = Zcapbanks_df(self.df_allcaps_selected)
        zcapbanks_obj2 = Zcapbanks(self.zcapbanks_obj.listofbanks)
        self.currents_capsfftobj = Current_capbanks_FFT(self.currents_buckfftobj,zcapbanks_obj2)
        self.i_zbank_ac_obj = Current_capbanks_waveforms(self.currents_capsfftobj)
        self.vcap_obj = Vripple_capbanks_waveforms(self.i_zbank_ac_obj,bw=10e6)

    def plots_multiple(self):
        fig, ax = plt.subplots(2,2)
        self.iobj.plot_waveforms(ax[0,0])
        self.i_zbank_ac_obj.plot_ibanks(ax[1,1])
        self.vcap_obj.plot_vripple(ax[0,1])
        self.zcapbanks_obj.plot_Zbanks(ax[1,0])

class InputCapbanks_a_vs_b:
    def __init__(self,ia,ib,db_0603,vpp_baseline=500):
        self.ia = ia; self.ib = ib
        self.db_0603 = db_0603
        self.vpp_baseline = vpp_baseline
        self.vpp_a = self.ia.vcap_obj.vpp; self.vpp_b = self.ib.vcap_obj.vpp
        self.df_a = self.ia.df_allcaps_selected
        self.df_b = self.ib.df_allcaps_selected
        add_sound_factor(self.df_a,self.db_0603,self.vpp_a/self.vpp_baseline)
        add_sound_factor(self.df_b,self.db_0603,self.vpp_b/self.vpp_baseline)
        self.db_total_a = self.db_total(self.df_a)
        self.db_total_b = self.db_total(self.df_b)
        
    def db_total(self,df_set):
        etotal = sum([db_to_energy_ratio(db_bank) for db_bank in df_set['sound factor(dB)'].values.tolist() if db_bank])
        return round(energy_ratio_to_db(etotal),1)
        
    def plots_a_vs_b(self):
        fig,ax = plt.subplots(2,3)
        fig.suptitle(f'Capgroup A ({self.db_total_a}dB) vs B ({self.db_total_b}dB)')
        fig.subplots_adjust(hspace=1)
        self.ia.i_zbank_ac_obj.plot_ibanks(ax[0,0])
        self.ia.vcap_obj.plot_vripple(ax[0,1])
        self.ia.zcapbanks_obj.plot_Zbanks(ax[0,2])
        self.ib.i_zbank_ac_obj.plot_ibanks(ax[1,0])
        self.ib.vcap_obj.plot_vripple(ax[1,1])
        self.ib.zcapbanks_obj.plot_Zbanks(ax[1,2])
        return fig

    def export_to_excel(self,df_a,df_b,fig):
        def write_plot_to_excel(fig,filename):
            buffer = BytesIO()
            fig.savefig(buffer, format='png', bbox_inches='tight')
            wb = load_workbook(filename)
            ws = wb['sheet1']#.active
            img = Image(buffer)
            ws.add_image(img,'D20')
            wb.save(filename)
        
        folder_path = r'sim_results'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        currentdate = datetime.datetime.today().strftime ('_%y%m%d_%H%M%S')
        filename = f'capbanks_a_vs_b'
        filename = filename+currentdate+'.xlsx' 
        outfile = os.path.join(folder_path, filename)
        params=self.ia.iobj.lobj.ckt['ip']
        df_params = pd.DataFrame.from_dict(params,orient='index')
        with pd.ExcelWriter(outfile) as writer:
            df_params.to_excel(writer,header=True, sheet_name='params')
            df_a.to_excel(writer,header=True, sheet_name="sheet1")
            df_b.to_excel(writer,header=True, startrow=7, sheet_name="sheet1")
            #writer.close()
        write_plot_to_excel(fig,outfile)
        print(f"File saved to: {outfile}")