

Rcu_in = {'1oz':10e-3,'2oz':5e-3}
Rcu_out = {'1oz':2e-3,'2oz':1e-3}

tcoeff = 3500e-6


def p_cu(i_inp,i_dc,oz_cu,tamb):
    tmult = tcoeff*(tamb-25)
    return(i_inp**2*Rcu_in[oz_cu]+i_dc**2*Rcu_out[oz_cu])*(1+tmult)

