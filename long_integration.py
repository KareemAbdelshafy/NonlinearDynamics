""" EXAMPLE: ABC Reaction  system found in [7]

    Drew LaMar, January 2006
"""
import PyDSTool as dst
from PyDSTool import *
import  math
import  numpy
from matplotlib import pyplot as plt

Q = 1000
w1 = 1.0
w2=2*w1

mu1 = w1/Q
mu2 = 2 *w1 /Q
G6 =  0.05  
phi = 0;

F1= 5e-05
F2= 0.0004
F1= 5e-05
F2= 0.00015
#F1= 1.5e-04
#F2= 0.001

OMG = w1

pars = { 'G6': G6, 'w1': w1, 'w2': w2, 'sig_W': 0., 'sig_w': 0., 'phi':phi, 'F1': F1, 'F2': F2, 'mu1': mu1, 'mu2': mu2, 'OMG': OMG}

print('F1 =', F1, 'F2 =', F2, 'w1 =', w1 , 'w2 =', w2)


icdict = {'x1': 0, 'x1D': 0, 'x2': 0 , 'x2D': 0}
# Set up model
x1Dstr = '-1*(mu1*x1D + w1^2*x1 + G6 *x1*x2 - F1*cos(OMG*t+phi) )'
x1str = 'x1D'
x2Dstr = '-1*(mu2*x2D + w2^2*x2 + 0.5*G6 *x1*x1 - F2*cos(2*OMG*t) )'
x2str = 'x2D'


DSargs = PyDSTool.args(name='Final1')
DSargs.pars = pars
DSargs.varspecs = {'x1D': x1Dstr, 'x1': x1str, 'x2D': x2Dstr, 'x2': x2str}
DSargs.ics = icdict
#DSargs.fnspecs  = {'r1': (['a1'], 'a1*a1') }
######
#DSargs.tdata = [0,25000]
DSargs.tdata = [0,50000]
#DSargs.algparams = {'init_step': 1}
DSargs.algparams = {'init_step': 2}
ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)   
traj = ode.compute('polarization')
pts  = traj.sample(dt=1)
x1= pts['x1'][20000:]
x2= pts['x2'][20000:]

plt.plot(pts['t'], pts['x1'])
plt.ylabel('x1')
plt.show()
plt.plot(pts['t'], pts['x2'])
plt.ylabel('x2')
plt.title(ode.name)                             # Figure title from model name
plt.show()
plt.plot(x1, x2,'.')
plt.xlabel(r'$a_1$', fontsize=20)
plt.ylabel(r'$a_2$', fontsize=20)
plt.show()
#####


