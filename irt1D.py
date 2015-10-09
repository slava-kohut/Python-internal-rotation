#!/usr/bin/python
from math import *
from numpy import linspace,zeros,arange,append,sort,amin,array
from scipy import sparse
from scipy.linalg import eigvalsh
from scipy.constants import N_A,hbar,R
from sys import argv,exit
from time import sleep
#
print 'irt1D (c) Sviataslau V. Kohut\n'
# grid size
k=float(argv[1]) 
ired=float(raw_input('Moment of inertia(*1e47 kg*m2): '))*10**(-47)
potential=raw_input('v(x): ')
print 'Wait for a while...'

code="""
def v(x):
    return %s
""" %potential
exec(code)
# elements
h=2*pi/k
A=-N_A*(hbar**2/(2*ired*(h**2)))
xi=linspace(-pi+h,pi,k)
vi=zeros(k)
for i in xrange(int(k)):
    vi[i]=(v(xi[i])-2*A)
# matrix
mat=zeros([k,k])
for i in xrange(int(k)):
 mat[i][i]=vi[i]
 if i-1<0:
  mat[i][i+k-1]=A
 else:
  mat[i][i-1]=A
 if i+1==k:
  mat[i][i-k+1]=A
 else:
  mat[i][i+1]=A
# diagonalization
levels=eigvalsh(mat)

filename=str(raw_input('Specify output file name: '))
outfile = open(filename, 'w')
title="""Moment of inertia: %.2fe-47 kg*m-2
Potential function: V(x)=%s\nLevels in J/mol:\n\n""" % (ired*10**47, potential)
outfile.write(title)

for k in xrange(255):
# only 255 levels will be present in an output file
    outfile.write('%8.2f' %levels[k])
    outfile.write('\n')
outfile.close()
# t/d functions
tdf=str(raw_input('Calculate contributions to thermodynamic functions? (y/n) '))
if tdf=='y':
     pass
else:
    exit(1)
filename=str(raw_input('Specify output file name: '))
outfile = open(filename, 'w')
title="""Moment of inertia: %.2fe-47 kg*m-2
Potential function: v(x)=%s\nContributions to thermodynamic functions:\n
    T       S      Cp     H/T   -G/T
    K     ----------J/K/mol----------\n""" % (ired*10**47, potential)
outfile.write(title)
temp=arange(50,1050,50);temp=append(temp,298.15);temp=append(temp,273.15);sort(temp)#array of temperatures
add=str(raw_input('Add extra values of temperature? default: [50:1000]+273.15,298.15 step:50 (y/n) '))
if add=='n':
     pass
else:
    while True: 
          value=str(raw_input('Input a value (c to continue): '))
          if value=='c':
              break
          else: 
              temp=append(temp,float(value))
ground=amin(levels)
for i in xrange(len(levels)):    
 # reducing levels to the ground state
 levels[i]=levels[i]-ground
# numerical derivatives
def Qnum(T):
    def qi(ei):
        return e**(-ei/(R*float(T)))
    qarray=zeros(len(levels))
    qarray=qi(levels)
    h=0.0005 #step
    q1=1./3.*sum(qarray)
    T=T-h
    qarray=zeros(len(levels))
    qarray=qi(levels)
    q0=1./3.*sum(qarray)
    T=T+2*h
    qarray=zeros(len(levels))
    qarray=qi(levels)
    q2=1./3.*sum(qarray)
    qln=(log(q2)-log(q0))/(2*h)
    qln2=(log(q2)-2*log(q1)+log(q0))/(h**2)
    return array([log(q1),qln,qln2]) #lnQ,d(lnQ)/dT,d2(lnQ)/dT2
# analytical derivatives
def Qan(T):
    def expqi(ei):
     return ei/(R*float(T)**2)
    def qi(ei):
     return e**(-ei/(R*float(T)))
    def dqidt(ei):
     return qi(ei)*(ei/(R*float(T)**2))
    def d2qidt2(ei):
     return dqidt(ei)*(expqi(ei)-2*T**-1)
    Q=zeros(len(levels))
    Q=sum(qi(levels))
    lnQdT=zeros(len(levels))  
    dQdT=sum(dqidt(levels))
    d2QdT2=zeros(len(levels))
    d2QdT2=sum(d2qidt2(levels))  
    return array([log(Q),dQdT/Q,d2QdT2/Q-dQdT**2/Q**2])         
   
def tdfunc(T,Q):
    S=R*Q(T)[0]+R*T*Q(T)[1]
    cp=2*R*T*Q(T)[1]+R*T**2*Q(T)[2]
    Hr=R*T*Q(T)[1]#H/T
    Gr=R*Q(T)[0]#-G/T
    return array([S,cp,Hr,Gr])

if argv[2]=='an':
 for T in sort(temp):
     outfile.write('%7.2f  %6.2f %6.2f %6.2f %6.2f' %(T,tdfunc(T,Qan)[0],tdfunc(T,Qan)[1],tdfunc(T,Qan)[2],tdfunc(T,Qan)[3]))
     outfile.write('\n')
 outfile.close()
elif argv[2]=='num':
 for T in sort(temp):
     outfile.write('%7.2f  %6.2f %6.2f %6.2f %6.2f' %(T,tdfunc(T,Qnum)[0],tdfunc(T,Qnum)[1],tdfunc(T,Qnum)[2],tdfunc(T,Qnum)[3]))
     outfile.write('\n')
 outfile.close()
else:
 print 'You entered illegal command-line argument!'
 sleep(3)
 exit(1)
     
     
     
    

    

    
    
    
