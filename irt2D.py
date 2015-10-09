from math import *
from numpy import empty,zeros,linspace,arange,append,amin,array
from math import pi
from scipy.constants import N_A,hbar,R
from numpy import zeros,sort
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix,csr_matrix
from time import sleep
from sys import argv

# command-line arguments
nx=int(argv[1])
ny=int(argv[2])
nlevels=int(argv[3])
vlevels=eval(argv[4])
#warning

if type(vlevels)!=str:
 print '4th command-line argument must be string, not %s' %(type(vlevels))
 sleep(3)
 exit(1)

# steps in x & y directions
hx=2*pi/nx
hy=2*pi/ny

Ired=float(raw_input('Moment of inertia(*1e47 kg*m2): '))*10**(-47)
Ax=-N_A*(hbar**2/(2*Ired*(hx**2)))
Ay=-N_A*(hbar**2/(2*Ired*(hy**2)))

xi=linspace(-pi+hx,pi,nx)
yk=linspace(-pi+hy,pi,ny)
# v(x,y)=v(x)+v(y)

potential=raw_input('V(x,y): ')
code="""def v(x,y):
  return %s""" %potential
exec(code)
 
# grid of potential values
grid_v=empty([nx,ny])
# calculating potential values at the gridpoints
for k in xrange(int(ny)):
 for i in xrange(int(nx)):
  grid_v[k][i]=v(xi[i],yk[k])

M=lil_matrix((int(nx*ny),int(nx*ny)))

#defining matrix elements...

print "forming the matrix..."

for k in xrange(int(ny)):
 for i in xrange(int(nx)):
  M[k*ny+i,k*ny+i]= grid_v[k][i]-2*(Ax+Ay) 
  if k-1<0:
   M[k*ny+i,(ny-1)*ny+i]= Ay
  else:
   M[k*ny+i,(k-1)*ny+i]= Ay
  if k+1==ny:
   M[k*ny+i,k*(ny-ny)+i]= Ay
  else:
   M[k*ny+i,(k+1)*ny+i]= Ay
  if i-1<0:
   M[k*ny+i,k*ny+(nx-1)]= Ax
  else:
   M[k*ny+i,k*ny+i-1]=Ax
  if i+1==nx:
   M[k*ny+i,k*ny+(nx-nx)]= Ax
  else:
   M[k*ny+i,k*ny+i+1]= Ax

M=csr_matrix(M)

print "done"
print """
diagonalization in process...
wait for a while..."""  
levels=sort(eigsh(M,nlevels,which=vlevels,tol=0.01,return_eigenvectors=False))

print 'energy levels calculated!i\n'

filename=str(raw_input('Specify output file name: '))
outfile = open(filename, 'w')
title="""A number of gridpoints: %d
Moment of inertia: %.2fe-47 kg*m-2
Potential function: V(x)=%s
Levels in J/mol:\n
-absolute--relative-\n""" % (nx*ny,Ired*10**47, potential)

outfile.write(title)
for k in xrange(len(levels)):
    outfile.write('%8.2f %8.2f' %(levels[k],levels[k]-amin(levels)))
    outfile.write('\n')
outfile.close()

tdf=str(raw_input('Calculate contributions to thermodynamic functions? (y/n) '))
if tdf=='y':
     pass
else:
    exit(1)
filename=str(raw_input('Specify output file name: '))
outfile = open(filename, 'w')
title="""Moment of inertia: ?
Potential function: v(x)=?\nContributions to thermodynamic functions:\n
    T       S      Cp     H/T   -G/T
    K     ----------J/K/mol----------\n"""
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
# reducing levels to the ground state
for i in xrange(len(levels)): 
    levels[i]=levels[i]-ground 

def Q(T):
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
    return array([log(Q),dQdT*Q**(-1),d2QdT2*Q**(-1)-dQdT**2*Q**(-2)])         
   
def tdfunc(T):
    S=R*Q(T)[0]+R*T*Q(T)[1]
    cp=2*R*T*Q(T)[1]+R*T**2*Q(T)[2]
    Hr=R*T*Q(T)[1]#H/T
    Gr=R*Q(T)[0]#-G/T
    return array([S,cp,Hr,Gr])

for T in sort(temp):
    outfile.write('%7.2f  %6.2f %6.2f %6.2f %6.2f' %(T,tdfunc(T)[0],tdfunc(T)[1],tdfunc(T)[2],tdfunc(T)[3]))
    outfile.write('\n')
outfile.close()
