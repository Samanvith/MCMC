import numpy as np
import math as m
import matplotlib.pyplot as plt
import pylab as pl 
from func import*
from cycler import cycle
#import Image
import matplotlib.ticker as mticker
#from constants import*


import math
from func import *


#file = np.loadtxt('jla_mub_0.txt')


#or om_ga in [0.3]:
#om_ga = 0.2
	#s = np.zeros((1,4))
	#z = np.linspace(0.010,2.00,1000)
	#a = A(z)
	#s = S(om_ga)
	#Et_a = eta(a,s)
	#del_f = delf(a,om_ga,z)
	#M_u = Mu(h,del_f)
	#plt.plot(z,M_u,'b-',label="")
	
	#print(s)



#print(z)
#print(s)

#print(del_f)

#print(Et_a)
#print(M_u)

#plt.plot(file[:, 0], file[:, 1], 'r-', label="observational datsampling a file")


n_sample = 10000

s_d = [0.01,0.01]

n_super = 31

param_val = np.empty([n_sample,3])

D_mu = np.empty(n_super)

param_val[0,:] = [np.random.uniform(), np.random.uniform(), 10**(-100)]

cov_mat = np.loadtxt('jla_mub_covmatrix.txt')

cov_mat = np.reshape(cov_mat,(n_super,n_super))

cov_Inv = np.linalg.inv(cov_mat)

Analytical_mU = np.empty(n_super)

Data = np.loadtxt('jla_mub_0.txt')

Z = Data[:,0]
m_u = Data[:,1]



def mu_analytical(om_ga,h,z):

	
	def A(z):
		return (1.0)/(1+z)
	#a = (1.0)/(1+z)


	def S(om_ga):
		return (np.power((1.0-om_ga)/om_ga,1.0/3.0))
	#s = ((1-om_ga)/om_ga)**(1.0/3.0)


	def eta(a,s):
		f_t = 2.0*np.power((np.power(s,3.0)+1.0),1.0/2.0)
		s_d = 1.0/np.power(a,4.0)
		t_d = 0.1540*(s/(np.power(a,3.0)))		
		f_h = 0.4304*((np.power(s,2.0))/(np.power(a,2.0)))
		f_v = 0.19097*((np.power(s,3.0))/(a))
		s_h = 0.066941 * np.power(s,4.0)

		part = np.power((s_d - t_d + f_h + f_v + s_h),(-1.0/8.0)) 

		e_ta = f_t*part
		return e_ta


	def delf(a,om_ga,z):
		return (c/Ho) * (1.0+z) * (eta(1.0,S(om_ga)) - eta(A(z),S(om_ga)))

	#S = S(om_ga)
		
	#del_f = c/Ho * (1+z) * ((eta(1,om_ga) - eta(a,om_ga)))
	
	def Mu(h,del_f):
		mu = 25.0-5.0*np.log10(h)+5.0*np.log10(del_f)
		return mu
	

	a = A(z)
	s = S(om_ga)
	Et_a = eta(a,s)
	del_f = delf(a,om_ga,z)
	M_u = Mu(h,del_f)

	return Mu(h,del_f)	

#for z in Z:
	#om_ga = 0.3
	#s = np.zeros((1,4))
	#z = np.linspace(0.010,2.00,1000)
#a = A(z)
#s = S(om_ga)
#Et_a = eta(a,s)
#del_f = delf(a,om_ga,z)

	
	#M_u = np.empty(n_super)
#M_u = Mu(h,del_f)
	
	#plt.plot(z,M_u,'b-',label="")
	
	#print(s)

pl.figure(figsize=(16.0,9.0))


def loglikehood(om_ga,h):
	if(om_ga<=0 or h<= 0):
		Log_Like = 10**(-100)
	else :
		for k in range (n_super):
			#print(Z)
			Analytical_mU[k] = mu_analytical(om_ga,h,Z[k])
			#print(len(Analytical_mU))
			D_mu[k] = m_u[k] - Analytical_mU[k]
			#print(D_mu)

		Log_Like = -0.5*np.dot(D_mu,np.dot(cov_Inv,D_mu))
		#print(Log_Like)
	return Log_Like

param_val[0,2] = loglikehood(param_val[0,0],param_val[0,1])

for i in range (1,n_sample):
	L_Like1 = param_val[i-1,2]
	om_ga_1 = np.random.normal(param_val[i-1,0],s_d[0])
	h_1 = np.random.normal(param_val[i-1,1],s_d[1])

	#print(i,len(Z))
	L_like2 = loglikehood(om_ga_1,h_1)
	

	if(L_like2 > L_Like1):
		param_val[i,0] = om_ga_1
		param_val[i,1] = h_1
		param_val[i,2] = L_like2
	else:
		A = np.random.uniform()
		if(L_like2-L_Like1>np.log(A)):
			param_val[i,0] = om_ga_1
			param_val[i,1] = h_1
			param_val[i,2] = L_like2
		else:
			param_val[i,0] = param_val[i-1,0]
			param_val[i,1] = param_val[i-1,1]
			param_val[i,2] = L_Like1



N2 = int(m.floor(n_sample/10))           




pl.scatter(param_val[N2:,0],param_val[N2:,1],c= -param_val[N2:,2])     
pl.xlim(0.15,0.45)
pl.ylim(0.65,0.75)
pl.xlabel("$\Omega m$-->")
pl.ylabel("h")
#pl.grid()
pl.legend(loc=4)
pl.savefig('Astrostats_4.png')

pl.show() 

pl.figure(figsize=(16.0,9.0))
pl.hist(param_val[N2:,0],bins=30)                              
pl.xlabel("$\Omega m$-->")
pl.legend(loc=4)
pl.savefig('Astrostats_5.png')

pl.show()











#plt.xlabel('z')
#plt.ylabel('Mu')
#plt.title('Z v/s Mu')
#plt.legend()
#plt.legend()
#plt.grid()
#plt.show()