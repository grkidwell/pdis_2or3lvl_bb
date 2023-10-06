
    
#to omit ac winding losses, set K1=0 in the Lparams when calling this function
#only need ac winding loss calculations for vishay inductors, whose model
#only works in CCM.  So will break this model for vishay inductors
def dcr_temp(i_rms_dcm,ckt_params:dict,Lparams:dict,p_core,tempco):
    cp = ckt_params.copy()
    lp = Lparams.copy()
    dcr = lp['DCR'] #at 25C
    Rth = lp['Rth']
    idc = cp['Idc']
    ipp = cp['deltaV']*cp['t_for_deltaV']/lp['Lout']/1e-6
    #tempco = 1/(234.45+25) #assume that the provided DCR spec is at 25C
    k1 = lp['K1']
    fs = cp['fs_phase']
    tamb = cp['Tamb']
    
    #dcr_term = dcr*Rth*(idc**2+k1*ipp**2*fs**0.5)

    
        
    dcr_term = dcr*Rth*i_rms_dcm**2
    return abs(((25*tempco-1)*dcr_term-tamb-p_core*Rth)/(dcr_term*tempco-1))