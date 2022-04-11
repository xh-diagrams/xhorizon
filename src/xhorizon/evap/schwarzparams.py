

import numpy as np
import matplotlib.pyplot as plt


"""
Consider an evaporating black hole with maximum mass M.
Suppose it evaporates by emitting Nevap shells of Hawking radiation.
This is described by a sequence of 
    Nevap+1 regions (first with m=M, last with m=0) with parameters m_j, du_j, dv_j, rp_j, rf_j
    Nevap shells (future boundary of each region) with parameters dmi, r_i, u_i, v_i
Suppose each is Schwarz so that
 	F(r) = 2m*( (r/2m-1) + log(r/2m-1) ).
Evap length is le.

j
M  ...      ... 0
0 | 1 | 2 | ...
  0   1   2  ...
  i

Constraints:
	R_j  = 2m_j
	dm_i = m_i+1 - m_i
	du_j = u_j - u_j-1
	dv_j = v_j - v_j-1

	dm_i < 0
	du_j = -3T/M m_j^2 dm_j-1 > 0
	dv_j - du_j = 2( F(rf_j) - F(rp_j) )

	rf_j = Rh(m_j)+le
	rp_j = Rh(m_j-1)+le

"""


def F(r,m):
	return r + 2.*m * np.log(np.abs(r/(2.*m) - 1))

def R(m):
	return 2.*m



def SSduvm(Nevap=1e4, Tevap=1e4, M=.5, le=1e-4):
	"""
	"""
	##
	eps = 1e-9
	s = np.linspace(eps,1.+eps,Nevap+1)[:-1]
	s = np.append(.5*eps,s)
	## u
	u  = Tevap * s
	## m(u)
	m = M * (1. - u/Tevap)**(1./3.)
	## r(u)
	r = 2.*m + le
	## v(u)
	v = u + 2.*F(r,m)
	## du dv dm
	du = u[1:]-u[:-1]
	dv = v[1:]-v[:-1]
	## mself mprev
	mself = m[1:]
	mprev = m[:-1]
	## drstar
	drstar = F(R(mself)+le,mself) - F(R(mprev)+le,mself)
	## zero
	zero = dv-du-2.*drstar
	## print 
	print("|dv-du-2.*drstar| (mean) = %s"%(np.mean(np.abs(zero))))
	print("|dv-du-2.*drstar| (total) = %s"%(Nevap*np.mean(np.abs(zero))))
	## return
	return mself, du, dv


#SSduvm()


