#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys
from qc import *
import tabulate

FILE_CPU = 'table_cpu'
PATH_TEX = '../doc/%s.tex'
PATH_CSV = '../doc/csv/%s.csv'

tabulate.LATEX_ESCAPE_RULES = {}

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
		return (mean, var, steps)

	def time_func(self, code, setup=''):
		exec(setup)
		tic = time.clock()
		exec(code)
		toc = time.clock()
		return toc-tic

class QProfilerCPU:
	def __init__(self):
		self.profile_cpu()

	def tic(self, d):
		d['tic'] = time.clock()

	def toc(self, d, e):
		d[e] = time.clock() - d['tic']
		d['tic'] = time.clock()


	def profile_cpu(self):
		Nmin = 2
		Nmax = 10
		Tmax = 60 #seconds
		Rmax = 100
		d = {}
		s = {}
		m = {}
		m2 = {}
		v = {}

		measures = ['qc0', 'qcf', 'm', 'cc', 'all']

		table = []
		for n in range(Nmin, Nmax+1):
			timming = time.clock()
			print(n)
			N = 2*n
			row = [n, N]
			for k in measures:
				s[k] = np.zeros(Rmax)

			for R in range(Rmax):
				config={"bits":n}

				self.tic(d)
				qs = QSimon(config)
				self.toc(d, measures[0])
				func = qs.random_config(n)
				st = qs.run(func)
				self.toc(d, measures[1])
				qm = QMeasure(st, np.arange(n, n*2))
				self.toc(d, measures[2])
				cp = CPostprocess(n, qm)
				cp.run()
				self.toc(d, measures[3])

				
				d['all'] = d['qc0'] + d['qcf'] + d['m'] + d['cc']

				for k in measures:
					s[k][R] = d[k]

				s['R'] = R+1
				if(timming + 60*5 < time.clock()): break

				#print(d)
			
			row.append(s['R'])

			for k in measures:
				array = s[k][0:R-1]
				m[k] = np.mean(array)
				m2[k] = np.log2(1+np.mean(array))/(2**N)
				v[k] = np.var(array)

			for h in measures:
				row.append(m[h])
				row.append(v[h])

			headers = ['n', 'N', 'R']
			for h in measures:
				headers.append(h + '-mean')
				headers.append(h + '-var')

			#print(m)
			#print(m2['all'])
		
			table.append(row)

		float_fmt = ".3e"
		table_fmt = 'simple'
		print(tabulate.tabulate(table, headers, floatfmt=float_fmt,
			tablefmt=table_fmt))

		#print(table)
		self.save_csv(table, headers, FILE_CPU)

	def save_csv(self, t, h, fn):
		str_header = ','.join(h)
		arrt = np.asarray(t)

		np.savetxt(PATH_CSV % fn, arrt, header=str_header, delimiter=",",
			comments='')

	def save_tex(self, t, h, fn):
		float_fmt = ".2f"
		table_fmt = 'latex_booktabs'
		str_table = tabulate.tabulate(t, h, floatfmt=float_fmt,
			tablefmt=table_fmt)

		self.save_file(str_table, PATH_TEX % fn)

	def save_file(self, s, fn):
		f = open(fn, 'w')
		f.write(s)
		f.close()


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

QProfilerCPU()
