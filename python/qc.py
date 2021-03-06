#!/usr/bin/env python3

# pyqsim - Quantum simulator of Simon algorithm.
# Copyright (C) 2016 Rodrigo Arias <rodrigo.arias@udc.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
		self.sz = {}
		self.size_all = -1
		self.U = None
		self.H = None
		self.phi1 = None
		self.qstate = None

		self.bits = self.config['bits']
		self.simon_init(self.bits)

	def simon_init(self, bits):
		x = ket(0, bits)
		y = ket(0, bits)

		H = hadamard(bits)
		I = identity(bits)
		#HI = tensor(H, I)
		#print("HI : {}".format(self.size_qobj(HI)))
		self.H = H

		phi0 = tensor(x, y)
		#phi1 = HI * phi0
		phi1 = self.HI(phi0)

		# Save some initial work
		self.phi1 = phi1
		self.qstate = None
		self.update_size()

	def size_qobj(self, o):
		if o == None: return 0
		data = o.data
		ptr = data.indptr.nbytes
		ind = data.indices.nbytes
		d = data.indices.nbytes
		return ptr + ind + d

	def update_size(self):
		self.sz = {}
		self.size_all = -1

		sz = self.sz
		sz['H'] = self.size_qobj(self.H)
		sz['U'] = self.size_qobj(self.U)
		sz['phi1'] = self.size_qobj(self.phi1)
		sz['qstate'] = self.size_qobj(self.qstate)
		sz['all'] = np.sum(list(sz.values()))
		#print('Updated size: {}'.format(str(sz)))

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
		phi1 = self.phi1
		f = params['f']

		self.U = self.build_U(f)
		#print('size U = %d' % self.size_qobj(U))
		phi2 = self.U * phi1
		#print(phi2.data.nnz, np.prod(phi2.data.shape))
		phi3 = self.HI(phi2)
		self.qstate = phi3

		self.update_size()

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
		self.sizes = {}
		self.prepare(qst, lines)

	def prepare(self, state, lines):
		n = len(lines)
		nz = state.data.nonzero()[0]
		ai = state.data[nz,0].T.toarray()[0]
		pi = np.array(np.real(np.power(ai, 2)))

		v = np.zeros(len(nz), dtype=int)
		for i in range(len(lines)):
			l = lines[i]
			v |= (nz & (1<<l)) >> (l-i)

		#print(v)
		#vp = np.zeros(len(pi))
		#for i in range(len(pi)):
		#	vp[i] = pi[v==i].sum()

		# O(2^(2^N)) -> O(2^N) !!!
		#vp = np.array([pi[v==i].sum() for i in range(len(pi))])
		vp = np.array([pi[v==i].sum() for i in range(2**n)])

		self.vp = vp[vp!=0]
		self.vn = np.arange(len(vp))[vp!=0]
		#print(self.vp)
		#self.update_size()

	def collapse(self):
		return int(np.random.choice(self.vn, p=self.vp))

	def update_size(self):
		self.sizes['vp'] = self.vp.nbytes
		self.sizes['vn'] = self.vn.nbytes
		#print(self.sizes)

class QProfiler:
	def __init__(self, steps=100, timeout=60*5):
		self.timeout = timeout
		self.steps = steps

	def start(self, code, setup='', title=None, fresh=True):
		timing = np.zeros(self.steps)
		tic = time.clock()
		if fresh:
			for i in range(self.steps):
				timing[i] = self.time_func(code, setup)
				toc = time.clock()
				if toc-tic > self.timeout: break
		else:
			exec(setup)
			for i in range(self.steps):
				ftic = time.clock()
				exec(code)
				ftoc = time.clock()
				timing[i] = ftoc-ftic
				toc = time.clock()
				if toc-tic > self.timeout: break

		timing = timing[0:i+1]
		mean = np.mean(timing)
		var = np.var(timing)
		steps = len(timing)
		title_str = ''
		if title: title_str = '%s: ' % title
		print("{}Mean {:.4e}, var {:.4e}, steps {}"
			.format(title_str, mean,var,steps))
		return (mean,var,steps)

	def time_func(self, code, setup=''):
		exec(setup)
		tic = time.clock()
		exec(code)
		toc = time.clock()
		return toc-tic

class CPostprocess():
	def __init__(self, bits, qm):
		self.bits = bits
		self.qm = qm

	def bitfield(self, a):
		bits = self.bits
		r = np.zeros(bits, dtype='int')
		for i in range(bits):
			if a & 1<<i:
				r[bits-i-1] = 1

		return r

	def solve(self, s):
		n = 2**self.bits
		for i in range(0, n):
			ok = True
			for j in s:
				x = i & j
				xbin = self.bitfield(x)
				if(list(xbin).count(1) % 2 == 1):
					ok = False
					break
			if not ok: continue
			#print('Solution {}'.format(i))

	def run(self):
		'Ejecuta la simulación de forma aleatoria'

		lin_ind = set()
		accum = set()
		steps = 0
		while True:
			if len(lin_ind) >= self.bits - 1: break
			y = self.qm.collapse()
			steps += 1
			if (y == 0) or (y in accum):
				#print('Vector miss {}'.format(y))
				continue
			tmp = accum.copy()
			for v in tmp:
				accum.add(v^y)
			accum.add(y)
			lin_ind.add(y)
			#print(lin_ind)

		#print('{} steps'.format(steps))
		#self.solve(lin_ind)
		return steps

	def profile(self, loops=100, timeout=60, title=None):
		vstep = np.zeros(loops)
		tic = time.clock()

		for i in range(loops):
			vstep[i] = self.run()
			toc = time.clock()
			if(toc - tic > timeout): break

		vstep = vstep[0:i+1]
		mean = np.mean(vstep)
		var = np.var(vstep)
		steps = len(vstep)

		title_str = ''
		if title: title_str = '%s: ' % title
		print("{}Mean {:.3e}, var {:.3e}, steps {}"
			.format(title_str, mean,var,steps))
		return (mean, var, steps)



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

	print("Profile of QMeasure.collapse()")
	for i in range(Nmin, Nmax+1):
		setup = ('qs = QSimon({"bits":%d});'%i) + \
				'config = qs.random_config(%d);'%i + \
				'st = qs.run(config);' + \
				'qm = QMeasure(st, np.arange(%d, %d*2))'%(i,i)
		qp.start('qm.collapse()', setup, title='N = %d'%i, fresh=False)

	print("Profile of np.random.choice()")
	for i in range(Nmin, 60):
		qp.start('np.random.choice(2**(%d-1))'%i, '', title='N = %d'%i)

def profile_mem():
	Nmin = 2
	Nmax = 10
	R = 10

	print("Profile of memory QSimon.init()")

	mems = list()

	for i in range(Nmin, Nmax+1):
		print("N = {}".format(i))
		mean_sizes = {}
		for r in range(R):
			N = i
			qs = QSimon({"bits":N})

			config = qs.random_config(N)
			#print(config)
			st = qs.run(config)
			sizes = qs.sz
			sizes['N'] = N
			mems.append(sizes)
			print(config)
			#print(sizes)
			#qm = QMeasure(st, np.arange(N, N*2))
			#cp = CPostprocess(N, qm)
			#cp.profile(10000, title='N = %d'%N)
			for k in sizes.keys():
				if not k in mean_sizes:
					mean_sizes[k] = sizes[k]
				else:
					mean_sizes[k] += sizes[k]
		print("N = {}, mean sizes: {}".format(i, mean_sizes))

	return mems


def main():
	#N = 10
	#if len(sys.argv) == 2: N = int(sys.argv[1])
	#print("N = %d" % N)

	#N = 2
	#qs = QSimon({"bits":2})
	#st = qs.run({'f':np.array([0, 0, 1, 1])})
	#qm = QMeasure(st, np.arange(N, N*2))
	#cp = CPostprocess(2, qm)
	#cp.run()

	#N = 9
	for N in range(2, 7):
		qs = QSimon({"bits":N})
		config = qs.random_config(N)
		#print(config)
		st = qs.run(config)
		qm = QMeasure(st, np.arange(N, N*2))
		cp = CPostprocess(N, qm)
		cp.run()
		cp.profile(10000, title='N = %d'%N)

	#st = qs.run({'f':np.array([1, 0, 1, 0])})
	##config = qs.random_config(N)
	#print(config)
	##st = qs.run(config)
	##qm = QMeasure(st, np.arange(N, N*2))
	#print(qm.vp, qm.vn)
	#rnd = np.random.choice(qm.vn, p=qm.vp, size=100)
	##print(qm.collapse())

#m = profile_mem()
#profile()
main()
