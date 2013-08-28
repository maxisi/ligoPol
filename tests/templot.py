import templates
import matplotlib.pyplot as plt
import numpy as np
import sidereal as sd

reload(templates)

days=1
t = np.array([x+ 630720013 for x in range(0, int(days*sd.ss), 60)])


det = 'H1'

kind = 'GR'
pdifs = ['0', 'p', 'm']

phi0=5.


plt.figure()

for p in pdifs:

    sig = templates.Signal('H1', 'J0534+2200', 'GR', p, t)
    
    s = sig.simulate(sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
    
    s.real.plot(label=p)
    
plt.title('Crab LHO GR signal (Re) $\phi_0=$' + str(phi0) )
plt.legend(numpoints=1, loc=4)


plt.figure()

for p in pdifs:

    sig = templates.Signal('H1', 'J0534+2200', 'GR', p, t)
    
    s = sig.simulate(sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
    
    s.imag.plot(label=p)
    
plt.title('Crab LHO GR signal (Im) $\phi_0=$' + str(phi0) )
plt.legend(numpoints=1, loc=4)


plt.figure()

for p in pdifs:

    sig = templates.Signal('H1', 'J0534+2200', 'GR', p, t)
    
    s = sig.simulate(sig.response.src.param['POL'], sig.response.src.param['INC'], phase=phi0)
    
    abs(s).plot(label=p)
    
plt.title('Crab LHO' + kind + ' signal amplitude $\phi_0=$' + str(phi0) )
plt.legend(numpoints=1, loc=4)

