#!/usr/bin/env python3

#from pyqtgraph.Qt import QtGui, QtCore
#import pyqtgraph as pg
import numpy as np
import qutip
import scipy.sparse.coo as coo

def alpha1(x):
	if x > 1: return 1
	return 0

def alpha2(x):
	t = [1,2,4,3]
	v = x % 4
	return t[v]

def alpha3(x):
	v = [4,1,5,7,1,4,7,5]
	return v[x]

def alpha4(x):
	v = [2,3,2,4,3,5,4,5,3,4,2,3,2,3,2,4,3,5,4,5,3,4,2,3,2,3,2,4,3,5,4,5,3,4,2,3]
	if x <= 35:
		return v[x]
	else:
		return 0

def build_alpha(bits, alpha):
	m = 2**bits
	x = np.zeros(m)
	alpha_x = np.zeros(m)
	for i in range(m):
		x[i] = i
		alpha_x[i] = alpha(i)
	return (x,alpha_x)

def build_pj(a, bits):
	n = a.size
	pj = np.zeros(n)
	for i in range(n):
		pj[i] = i << bits | int(a[i])
	return pj

def build_pi(a, bits):
	n = a.size
	pi = np.zeros(n)
	for i in range(n):
		pi[i] = i << bits
	return pi

def build_vus(pi, pj):
	n = pi.size
	m = n**2
	us = np.zeros(m)
	for i in range(m):
		us[i] = i
	for i in range(n):
		us[int(pi[i])] = pj[i]
		us[int(pj[i])] = pi[i]
	return us

def hadamard(q):
	n = q.size
	bits = int(np.log2(n))

def collapse(a, b):
	v = np.random.random_sample()
	if v < a: r = 0
	else: r = 1

	print("collapse ({} {})T = {}".format(a, b, r))
	
	return r

def measure(c, k, n):
	p_alpha = 0
	p_beta = 0
	for q in range(2**n):
		if (int(q) & int(1<<k)) == 0:
#			print("alpha {}".format(q))
			p_alpha += np.absolute(c[q][0][0])**2
		else:
#			print("beta {}".format(q))
			p_beta += np.absolute(c[q][0][0])**2
	alpha = np.sqrt(p_alpha)
	beta  = np.sqrt(p_beta)
#	print("alpha {}".format(alpha))
#	print("beta {}".format(beta))

	x = np.zeros(2**n, dtype=complex)
	for q in range(2**n):
		if (int(q) & int(1<<k)) == 0:
			if alpha != 0:
				x[q] = c[q][0][0]/alpha
		else:
			if beta != 0:
				x[q] = c[q][0][0]/beta

	r = collapse(p_alpha, p_beta)
#	print("Measure {} = {}".format(k, r))

	for q in range(2**n):
		if q&(2**k) != 0: w = 1
		else: w = 0
		if w != r:
			x[q] = 0
	
	st = qutip.Qobj(x)

	return (st, r)

def measures(c, v, n):
	st = c
	result = []
	for i in range(len(v)):
		(st, r) = measure(st, v[i], n)
#		print(st)
		result.append(r)

	return (st, result)

def hadamard(bits):
	return qutip.Qobj(qutip.hadamard_transform(bits).data)

def identity(bits):
	return qutip.Qobj(qutip.identity(2**bits).data)

def h_i(bits):
	h = hadamard(bits)
	i = identity(bits)
	return qutip.Qobj(qutip.tensor(h, i).data)

def build_us(bits, fun):
	n = 2**(bits*2)
	(x, alpha_x) = build_alpha(bits, fun)
	pi = build_pi(alpha_x, bits)
	pj = build_pj(alpha_x, bits)
	us = build_vus(pi, pj)

	max_n = len(us)
	data = [1]*max_n
	vi = us
	vj = np.arange(0, max_n)
	Us_mat = coo.coo_matrix((data, (vi, vj)), shape=[n,n])
	Us = qutip.Qobj(Us_mat)
	return Us

def simon(bits, fun):
	n = 2**(bits*2)

	qx = qutip.ket([0]*bits)
	qy = qutip.ket([0]*bits)

	phi0 = qutip.Qobj(qutip.tensor(qx, qy).data)
	phi1 = h_i(bits) * phi0
	phi2 = build_us(bits, fun) * phi1
	phi3 = h_i(bits) * phi2
	measure_i = np.arange(bits, 2*bits)
#	print(phi3)
#	print([(i,e) for i, e in enumerate(phi3) if e != 0])
	(phi4, r) = measures(phi3, measure_i, 2*bits)
#	print(phi3)
	return r

def shifting(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out


f = [alpha1, alpha2, alpha3, alpha4]
b = [2, 4, 3, 5]
fi = 3
bits = b[fi]
fun = f[fi]
N = 500

count = np.zeros(2**bits)
med = np.zeros(bits)
np.set_printoptions(precision=3)
#app = QtGui.QApplication([])
#win = pg.GraphicsWindow()
#pl = win.addPlot()
x = np.arange(2**bits)
for i in range(N):
	r = simon(bits, fun)
	med += r
	v = shifting(r)
	count[v] += 1
	mp = med / (i+1)
	print("measured {}, count {}".format(v, count))
#	for j in range(len(count)):
#		pl.setData(j, count[j])

print(count)

