
Rcu_in = 0.5e-3
Rcu_out = 1e-3

tcoeff = 3500e-6


def p_cu(i_inp,i_dc,tamb):
    tmult = tcoeff*(tamb-25)
    return(i_inp**2*Rcu_in+i_dc**2*Rcu_out)*(1+tmult)