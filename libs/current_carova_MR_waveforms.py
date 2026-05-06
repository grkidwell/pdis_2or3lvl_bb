import numpy as np
import scipy as sp
import scipy.signal

class Current_carova_MR_waveforms:
    def __init__(self, ind_obj):
        self.lobj = ind_obj
        self.fs = self.lobj.ckt['ip']['fs']
        self.ts = 1/self.fs
        self.tarray = np.linspace(0,self.ts/2,512)
        #self.time_datapoints is created in 2cycle_ccm function below
        self.ind_wf_dcmpulse = self.create_ind_waveform_dcmpulse(self.tarray)
        self.i_ind = self.create_ind_waveform_2cycle_ccm()
        self.i_hs = self.create_hs_waveform()
        self.i_zin_ac = self.i_hs-np.mean(self.i_hs)
        
    def create_ind_waveform_dcmpulse(self,t):
        tstate = self.lobj.ckt['t_state']; deltaV = self.lobj.ckt['deltaV']; lout = self.lobj.ckt['ip']['lout']['value(uH)']*1e-6
        tstart = {1:0, 2:tstate[1], 3:tstate[1]+tstate[2]}      
        tend = {state:tstart[state]+tstate[state] for state in tstate.keys()}
        ipp = {state:round(dv*tstate[state]/lout,3) for state,dv in deltaV.items()}   
        i_start = {1:0, 2:ipp[1], 3:ipp[1]+ipp[2]}
        def waveform(t):
            def conditions(t):
                return [t < tend[1], (t>=tend[1]) & (t<tend[2]), t>=tend[2]]           
            def functions():
                return [lambda x:i_start[1]+deltaV[1]*(x-tstart[1])/lout, 
                        lambda x:i_start[2]+deltaV[2]*(x-tstart[2])/lout, 
                        lambda x:i_start[3]+deltaV[3]*(x-tstart[3])/lout]
            return np.piecewise(t,conditions(t),functions())
        return waveform(t)   
    def create_ind_waveform_2cycle_ccm(self):
        t_2ndperiod = self.tarray+self.ts/2
        self.time_datapoints = np.concatenate((self.tarray,t_2ndperiod))
        wf_2cycles = np.concatenate((self.ind_wf_dcmpulse,self.ind_wf_dcmpulse))
        avg = float(round(self.ind_wf_dcmpulse.mean(),3))
        return wf_2cycles-avg+self.lobj.ckt['Idc']

    def create_hs_waveform(self):
        tstate = self.lobj.ckt['t_state']
        d_hs_1stcycle = tstate[1]+tstate[2]
        D_pulse_1stcycle = sp.signal.square(2*np.pi*self.tarray,d_hs_1stcycle)/2+0.5
        d_hs_2ndcycle = tstate[1]
        D_pulse_2ndcycle = sp.signal.square(2*np.pi*self.tarray,d_hs_2ndcycle)/2+0.5
        D_pulse_2cycles = np.concatenate((D_pulse_1stcycle,D_pulse_2ndcycle))
        return D_pulse_2cycles*self.i_ind

    def plot_waveforms(self,ax):
        time_datapoints = self.time_datapoints
        ax.plot(time_datapoints,self.i_ind,label='inductor')
        ax.plot(time_datapoints,self.i_zin_ac,label='input cap total')
        ax.legend(loc='lower left',bbox_to_anchor=(0.0,1.02));ax.grid(True,color='lime',linewidth=.2)
    
