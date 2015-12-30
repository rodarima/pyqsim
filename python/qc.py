#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys

def ket(st, bits):
	return qutip.basis(2**bits, st)

def hadamard(bits):
	return qutip.Qobj(qutip.hadamard_transform(bits).data)

def identity(bits):
	return qutip.identity(2**bits)

def tensor(a, b):
	return qutip.Qobj(qutip.tensor(a, b).data)

class QCircuit:
	def __init__(self, c):
		self.config = c
		self.init()

	def init(self): pass

	def run(self, params): pass


class QSimon(QCircuit):
	def init(self):
		self.bits = self.config['bits']
		self.simon_init(self.bits)

	def simon_init(self, bits):
		x = ket(0, bits)
		y = ket(0, bits)

		H = hadamard(bits)
		I = identity(bits)
		#HI = tensor(H, I)
		self.H = H

		phi0 = tensor(x, y)
		#phi1 = HI * phi0
		phi1 = self.HI(phi0)

		# Save some initial work
		self.phi1 = phi1
		self.qstate = phi1

	def HI(self, state):
		bits = self.bits

		nz = state.data.nonzero()[0]
		ai = state.data[nz,0].T.toarray()[0]
		iH = nz >> bits
		iL = nz % 2**bits

		tmp = qutip.zero_ket(2**(2*bits))
		for i in range(len(nz)):
			a = ai[i]
			ketH = self.H * ket(iH[i], bits)
			ketL = ket(iL[i], bits)
			tmp += a * tensor(ketH, ketL)

		return tmp

	def run(self, params):
		# Start from the precomputed state
		phi1 = self.qstate
		f = params['f']

		U = self.build_U(f)
		phi2 = U * phi1
		#print(phi2.data.nnz, np.prod(phi2.data.shape))
		phi3 = self.HI(phi2)
		self.qstate = phi3
		return phi3

	def build_U(self, f):
		bits = self.bits
		n = 2**(2*bits)
		if(bits>32):
			print("bits = {}, exiting".format(bits))
			exit(1)
		b = np.arange(n)
		bH = b >> bits
		bL = b % 2**bits
		aL = f[bH] ^ bL
		a = bH << bits | aL
		data = [1]*n
		U = coo.coo_matrix((data, (a, b)), shape=[n,n])
		return qutip.Qobj(U)

	def random_config(self, bits):
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
		return {'f':alpha, 'period':s}

class QMeasure:
	def __init__(self, qst, lines):
		self.qst = qst
		self.prepare(qst, lines)

	def prepare(self, state, lines):
		nz = state.data.nonzero()[0]
		ai = state.data[nz,0].T.toarray()[0]
		pi = np.array(np.real(np.power(ai, 2)))

		v = np.zeros(len(nz), dtype=int)
		for i in range(len(lines)):
			l = lines[i]
			v |= (nz & (1<<l)) >> (l-i)

		vp = np.zeros(len(pi))
		for i in range(len(pi)):
			vp[i] = pi[v==i].sum()

		self.vp = vp[vp!=0]
		self.vn = np.arange(len(vp))[vp!=0]

	def collapse(self):
		return np.random.choice(self.vn, p=self.vp)

class QProfiler:
	def __init__(self, steps=100, timeout=60*5):
		self.timeout = timeout
		self.steps = steps

	def start(self, code, setup='', title=None):
		timing = np.zeros(self.steps)
		tic = time.clock()
		for i in range(self.steps):
			timing[i] = self.time_func(code, setup)
			toc = time.clock()
			if toc-tic > self.timeout: break

		timing = timing[0:i+1]
		mean = np.mean(timing)
		var = np.var(timing)
		steps = len(timing)
		title_str = ''
		if title: title_str = '%s: ' % title
		print("{}Mean {:.3e}, var {:.3e}, steps {}"
			.format(title_str, mean,var,steps))

	def time_func(self, code, setup=''):
		exec(setup)
		tic = time.clock()
		exec(code)
		toc = time.clock()
		return toc-tic

def profile():
	Nmin = 2
	Nmax = 9
	Tmax = 60 #seconds

	qp = QProfiler(timeout=Tmax)

	print("Profile of QSimon.init()")
	for i in range(Nmin, Nmax+1):
		qp.start('QSimon(config)',
			'config={"bits":%d}'%i,
			title='N = %d'%i)

	print("Profile of QSimon.run()")
	for i in range(Nmin, Nmax+1):
		setup = ('qs = QSimon({"bits":%d});'%i) + \
				'config = qs.random_config(%d)'%i
		qp.start('qs.run(config)', setup, title='N = %d'%i)

	print("Profile of QMeasure()")
	for i in range(Nmin, Nmax+1):
		setup = ('qs = QSimon({"bits":%d});'%i) + \
				'config = qs.random_config(%d);'%i + \
				'st = qs.run(config)'
		qp.start('QMeasure(st, np.arange(%d, %d*2))'%(i,i), setup, title='N = %d'%i)


def main():
	#N = 10
	#if len(sys.argv) == 2: N = int(sys.argv[1])
	#print("N = %d" % N)

	N = 2
	qs = QSimon({"bits":2})
	st = qs.run({'f':np.array([0, 0, 1, 1])})
	qm = QMeasure(st, np.arange(N, N*2))


	#st = qs.run({'f':np.array([1, 0, 1, 0])})
	##config = qs.random_config(N)
	#print(config)
	##st = qs.run(config)
	##qm = QMeasure(st, np.arange(N, N*2))
	#print(qm.vp, qm.vn)
	#rnd = np.random.choice(qm.vn, p=qm.vp, size=100)
	##print(qm.collapse())

#main()
profile()
