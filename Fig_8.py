""" EXAMPLE: ABC Reaction  system found in [7]

    Drew LaMar, January 2006
"""
import PyDSTool as dst
from PyDSTool import *
import  math
import  numpy
from matplotlib import pyplot as plt

Q = 1000
G1 = 1 #0.4247
G2 = 1 #1.6623
G3 = 1 #-192
w1 = math.sqrt(G3/G1)
w2=2*w1
G4 = 4 * G3 * G2 / G1

mu1 = w1/Q
mu2 = 2 *w1 /Q
G6 =  0.05 #0.05 #115

VFDC = 1
VFAC = 0.04 #0.02 #0.1 
VPDC = 0.1
VPAC = 5
F1= 5e-05
F2= 0.0005



pars = {'G1': G1, 'G2': G2, 'G3': G3, 'G4': G4, 'G6': G6, 'w1': w1, 'w2': w2, 'sig_W': 0., 'sig_w': 0., 'phi':-pi/2, 'F1': F1, 'F2': F2, 'mu1': mu1, 'mu2': mu2}

print('F1 =', F1, 'F2 =', F2, 'w1 =', w1 , 'w2 =', w2)


icdict = {'a1': -0.05, 'a2': 0.00005, 'b1':  0.00006 , 'b2': -0.025 }
# Set up model
b1str = '(2*w1*sig_W*a1 +w1*mu1*b1 +0.25*9*G6*G6/8/w2^2*a1*(a1*a1+b1*b1) +0.25*G6*G6/8/w1^2*a1*(a2*a2+b2*b2) -0.5*G6*(a1*a2+b1*b2) +F1*cos(phi))/(-2*w1)'
a1str = '(2*w1*sig_W*b1 -w1*mu1*a1 +0.25*9*G6*G6/8/w2^2*b1*(a1*a1+b1*b1) +0.25*G6*G6/8/w1^2*b1*(a2*a2+b2*b2) -0.5*G6*(a1*b2-b1*a2) +F1*sin(phi))/(2*w1)'
b2str = '(2*w2*(2*sig_W-sig_w)*a2 +w2*mu2*b2 +0.25*G6*G6/8/w1^2*a2*(a1*a1+b1*b1) -0.25*G6*(a1*a1-b1*b1) +F2)/(-2*w2)'
a2str = '(2*w2*(2*sig_W-sig_w)*b2 -w2*mu2*a2 + 0.25*G6*G6/8/w1^2*b2*(a1*a1+b1*b1) -0.5*G6*a1*b1)/(2*w2)'


DSargs = PyDSTool.args(name='Final1')
DSargs.pars = pars
DSargs.varspecs = {'b1': b1str, 'a1': a1str, 'b2': b2str, 'a2': a2str}
DSargs.ics = icdict
#DSargs.fnspecs  = {'r1': (['a1'], 'a1*a1') }
######
'''
DSargs.tdata = [0,25000]
DSargs.algparams = {'init_step': 1e-1}
ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)   
traj = ode.compute('polarization')
pts  = traj.sample(dt=100)

print(type(pts['a2']))
A22 = (pts['a2'])*(pts['a2']) + (pts['b2'])*(pts['b2'])
A2 = numpy.sqrt(A22[:])

A12 = (pts['a1'])*(pts['a1']) + (pts['b1'])*(pts['b1'])
A1 = numpy.sqrt(A12[:])

plt.plot(pts['t'], A2)
plt.ylabel('A2')
plt.show()
plt.plot(pts['t'], A1)
plt.ylabel('A1')
plt.show()

plt.plot(pts['t'], pts['a2'])
plt.ylabel('a2')
plt.show()
plt.plot(pts['t'], pts['b2'])
plt.ylabel('b2')
plt.show()
plt.plot(pts['t'], pts['a1'])
plt.ylabel('a1')
plt.show()
plt.plot(pts['t'], pts['b1'])
plt.ylabel('b1')                           # ...
#plt.ylim([0,65])                                # Range of the y axis
plt.title(ode.name)                             # Figure title from model name
plt.show()
#####

'''
testDS = Generator.Vode_ODEsystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['sig_W']
PCargs.StepSize = 80e-5
PCargs.MaxNumPoints = 400  #500
PCargs.MaxStepSize = 100e-5
PCargs.MinStepSize = 1e-6
PCargs.NumSPOut = 30;
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PCargs.auxpars = {'r1'}
#PCargs.fnspecs  = {'r1': (['a1'], 'a1*a1') }
#PCargs.MaxTestIters = 10
#PCargs.FuncTol = 1e-10
#PCargs.VarTol = 1e-8
#PCargs.TestTol = 1e-7

PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()
print('done in %.3f seconds!' % (clock()-start))

sol = PyCont['EQ1'].sol


A11 = numpy.sqrt(sol['a1'] *sol['a1'] + sol['b1']* sol['b1'])
A21 = numpy.sqrt(sol['a2'] *sol['a2'] + sol['b2']* sol['b2'])
sig_W1= sol['sig_W']


PyCont.display(('sig_W','a2'), stability=True)
PyCont.display(('sig_W','b2'), stability=True)
PyCont.display(('sig_W','a1'), stability=True)
PyCont.display(('sig_W','b1'), stability=True)
show()
plt.plot(sig_W1, A21)
plt.plot(sig_W1, A11)
show()
'''
'''
#######################################
#F2= 0

pars['F2'] = F2
icdict = {'a1': 0.0704, 'a2': -0.00005, 'b1':  -0.00006*2 , 'b2': 0.0016*2 }
DSargs = PyDSTool.args(name='Final1')
DSargs.pars = pars
DSargs.varspecs = {'b1': b1str, 'a1': a1str, 'b2': b2str, 'a2': a2str}
DSargs.ics = icdict
testDS = Generator.Vode_ODEsystem(DSargs)
PyCont = ContClass(testDS)
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['sig_W']
PCargs.StepSize = 80e-5
PCargs.MaxNumPoints = 400  #500
PCargs.MaxStepSize = 100e-5
PCargs.MinStepSize = 1e-6
PCargs.NumSPOut = 30;
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()
print('done in %.3f seconds!' % (clock()-start))

sol = PyCont['EQ1'].sol
#print(sol.info(2))

A12 = numpy.sqrt(sol['a1'] *sol['a1'] + sol['b1']* sol['b1'])
A22 = numpy.sqrt(sol['a2'] *sol['a2'] + sol['b2']* sol['b2'])
sig_W2= sol['sig_W']

############################################

pars['G6'] = 0
icdict = {'a1': 0.0704, 'a2': -0.00005, 'b1':  -0.00006*2 , 'b2': 0.0016*2 }
DSargs = PyDSTool.args(name='Final1')
DSargs.pars = pars
DSargs.varspecs = {'b1': b1str, 'a1': a1str, 'b2': b2str, 'a2': a2str}
DSargs.ics = icdict
testDS = Generator.Vode_ODEsystem(DSargs)
PyCont = ContClass(testDS)
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['sig_W']
PCargs.StepSize = 80e-5
PCargs.MaxNumPoints = 400  #500
PCargs.MaxStepSize = 100e-5
PCargs.MinStepSize = 1e-6
PCargs.NumSPOut = 30;
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PCargs.SaveEigen = True
PCargs.SaveJacobian = True
PyCont.newCurve(PCargs)

print('Computing curve...')
start = clock()
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()
print('done in %.3f seconds!' % (clock()-start))

sol = PyCont['EQ1'].sol
#print(sol.info(2))

A13 = numpy.sqrt(sol['a1'] *sol['a1'] + sol['b1']* sol['b1'])
A23 = numpy.sqrt(sol['a2'] *sol['a2'] + sol['b2']* sol['b2'])
sig_W3= sol['sig_W']

############################################


PyCont.display(('sig_W','a2'), stability=True)
PyCont.display(('sig_W','b2'), stability=True)
PyCont.display(('sig_W','a1'), stability=True)
PyCont.display(('sig_W','b1'), stability=True)
show()
'''
'''
a11, =plt.plot(sig_W1, A11, linewidth=2.0 , color = 'blue')
a12, = plt.plot(sig_W2, A12, linewidth=2.0, color = 'blue')
a13, = plt.plot(sig_W3, A13, linewidth=2.0, color = 'red')
plt.legend([a11,a13] , [r'$\alpha=0.05$',r'$\alpha=0$'], fontsize=20, loc= 'best')
plt.xlabel(r'$\sigma_\Omega$', fontsize=20)
plt.ylabel(r'$a_1$', fontsize=20)
plt.axis([-0.004, 0.004, 0, 0.21])
show()

sig_W3 = sig_W3[::20]
A23 = A23[::20]

a21, =plt.plot(sig_W1, A21, linewidth=2.0 , color = 'blue')
a22, = plt.plot(sig_W2, A22, linewidth=2.0, color = 'blue')
a23, = plt.plot(sig_W3, A23,'o', linewidth=2.0, color = 'red')
plt.xlabel(r'$\sigma_\Omega$', fontsize=20)
plt.legend([a21,a23] , [r'$\alpha=0.05$',r'$\alpha=0$'], fontsize=20, loc= 'best')
plt.ylabel(r'$a_2$', fontsize=20)
plt.axis([-0.004, 0.004, 0, 0.14])
show()
