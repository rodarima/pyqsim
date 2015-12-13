#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time

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
	if x <= 15:
		return v[x]
	else:
		return 0

def alpha5(x):
	v = [2, 3]
	return v[x]

def alpha6(x):
	v = [2, 3, 2, 4]
	return v[x]

def alpha7(x):
	v = [2, 3, 2, 4, 3, 5, 4, 5]
	return v[x]

def alpha632(x):
	v = [0,4,4,0,2,6,6,2]
	return v[x]

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

#	print("collapse ({} {})T = {}".format(a, b, r))
	
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
	bits = len(v)
	result = np.zeros(bits, dtype='int')
	for i in range(bits):
		(st, r) = measure(st, v[i], n)
#		print(st)
		result[bits-i-1] = r

	return (st, result)

def hadamard(bits):
	return qutip.Qobj(qutip.hadamard_transform(bits).data)

def identity(bits):
	return qutip.Qobj(qutip.identity(2**bits).data)

def h_i(bits):
	h = hadamard(bits)
	i = identity(bits)
	return qutip.Qobj(qutip.tensor(h, i).data)

def build_us(bits, alpha):
	n = 2**(bits*2)
	pi = build_pi(alpha, bits)
	pj = build_pj(alpha, bits)
	us = build_vus(pi, pj)

	max_n = len(us)
	data = [1]*max_n
	vi = us
	vj = np.arange(0, max_n)
	Us_mat = coo.coo_matrix((data, (vi, vj)), shape=[n,n])
	Us = qutip.Qobj(Us_mat)
	return Us

def check(bits, alpha, s):
	n = 2**bits
	for i in range(n):
		mirror = i ^ s
		if alpha[i] != alpha[mirror]:
			print("Error, la función no cumple las condiciones")
			print("alpha({}) = {} != alpha({}) = {}".format(
				i, alpha[i], mirror, alpha[mirror]))
			exit()

def simon(bits, alpha, s):
	n = 2**(bits*2)
	check(bits, alpha, s)

	qx = qutip.ket([0]*bits)
	qy = qutip.ket([0]*bits)

	phi0 = qutip.Qobj(qutip.tensor(qx, qy).data)
	phi1 = h_i(bits) * phi0
	#tic = time.clock()
	phi2 = build_us(bits, alpha) * phi1
	#print("Time build_us = {}".format(time.clock() - tic))
	phi3 = h_i(bits) * phi2
	measure_i = np.arange(bits, 2*bits)
#	print(phi3)
#	print([(i,e) for i, e in enumerate(phi3) if e != 0])
	tic = time.clock()
	(phi4, r) = measures(phi3, measure_i, 2*bits)
	print("Time measure = {}".format(time.clock() - tic))
#	print(phi3)
	return r

def shifting(bitlist):
	out = 0
	for bit in bitlist:
		out = (out << 1) | bit
	return out

def bitfield(a, bits):
	r = np.zeros(bits, dtype='int')
	for i in range(bits):
		if a & 1<<i:
			r[bits-i-1] = 1

	return r

def check_solution(eq, sol, bits):
	eqv = bitfield(eq, bits)
	solv = bitfield(sol, bits)
	n = len(eqv)
	r = np.zeros(n, dtype='int')
	for i in range(len(eqv)):
		r[i] = eqv[i]*solv[i]

def solve(v, bits):
	n = 2**bits
	s = np.zeros(n)
	for i in range(0, n):
		for j in range(n):
			if v[j] > 0:
				x = i&j
				xbin = bitfield(x, bits)
				ones = list(xbin).count(1) % 2
#				print("{} and {} = {} -> {}".format(i, j, xbin, ones))
				s[i]+=ones
	sol = []
	for i in range(n):
		if s[i] == 0:
			sol.append(i)

#	print("Solutions {}".format(sol))
	return sol

def random_alpha(bits):
	n = 2**bits
	alpha = np.zeros(n, dtype='int')
	alpha.fill(-1)
	s = np.random.randint(1, n)
	for i in range(n):
		if(alpha[i] == -1):
			r = np.random.randint(0, n)
			while r in list(alpha):
				r = np.random.randint(0, n)
			alpha[i] = r
			alpha[i^s] = r
	return (alpha, s)

f = [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha632]
s = [1, 4, 5, 12, 1, 2, 12, 3]
b = [2, 4, 3, 4, 1, 2, 3, 3]
fi = 7

def test_simon(bits, max_N):

	(alpha, s) = random_alpha(bits)

	count = np.zeros(2**bits, dtype='int')
	med = np.zeros(bits)
	np.set_printoptions(precision=3)
	x = np.arange(2**bits)

	i = 0
	while True:
#		print("--- Iteration {} ---".format(i))
		#tic = time.clock()
		r = simon(bits, alpha, s)
		#t_simon = time.clock() - tic
		med += r
		v = shifting(r)
		count[v] += 1
		#tic = time.clock()
		sol = solve(count, bits)
		#t_solve = time.clock() - tic
		mp = med / (i+1)
		#print("Time simon = {}, time solve = {}".format(t_simon, t_solve))
#		print("measured {}, count {}".format(v, count))
#		print("Solutions {}".format(sol[1:]))
		if len(sol) == 2: break
		i+=1


	if len(sol) == 2:
#		print("Period is {} = {}".format(sol[1], bitfield(sol[1], bits)))
#		print("Expected  {} = {}".format(s, bitfield(s, bits)))
		if s!=sol[1]:
			print("ERROR period incorrect!")
			print((alpha, s))
			exit()
	else:
		print("ERROR, period not found")
		print((alpha, s))
		exit()

#	print("Solution found after {} steps".format(i+1))
	return i+1


def main():
	bits = 4
	N = 10
	max_it = 1000
	steps = 0
	step_v = np.zeros(N)
	step_t = np.zeros(N)

	print("Simulation of Simon quantum algorithm.")
	print("bits = {}, simulations = {}".format(bits, N))
	print("N\tsteps\ttime\ttime/steps")

	for i in range(N):
		tic = time.clock()
		p = test_simon(bits, max_it)
		toc = time.clock()
		t = toc - tic
		step_v[i] = p
		step_t[i] = t

		print("{}\t{}\t{:.3f}\t{:.3f}".format(i, p, t, t/p))

	mean_n = np.mean(step_v)
	var_n = np.var(step_v)
	mean_t = np.mean(step_t / step_v)
	var_t = np.var(step_t / step_v)

	print("Analysis of simulated data.")
	print("Mean steps = {}".format(mean_n))
	print("Var steps = {}".format(var_n))
	print("Mean time/step = {}".format(mean_t))
	print("Var time/step = {}".format(var_t))

main()
