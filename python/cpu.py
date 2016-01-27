#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys
from qc import *
import tabulate

FILE_QC0 = 'cpu_qc0'
FILE_QCF = 'cpu_qcf'
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

def tic(d):
	d['tic'] = time.clock()

def toc(d, e):
	d[e] = time.clock() - d['tic']
	d['tic'] = time.clock()


def profile_cpu():
	Nmin = 2
	Nmax = 6
	Tmax = 60 #seconds
	Rmax = 10
	d = {}
	s = {}
	m = {}
	m2 = {}
	v = {}

	measures = ['qc0', 'qcf', 'm', 'cc', 'all']
	headers = ['N'] + measures

	table = []
	for N in range(Nmin, Nmax+1):
		row = [N]
		for k in measures:
			s[k] = np.zeros(Rmax)

		for R in range(Rmax):
			config={"bits":N}

			tic(d)
			qs = QSimon(config)
			toc(d, measures[0])
			func = qs.random_config(N)
			st = qs.run(func)
			toc(d, measures[1])
			qm = QMeasure(st, np.arange(N, N*2))
			toc(d, measures[2])
			cp = CPostprocess(N, qm)
			cp.run()
			toc(d, measures[3])

			
			d['all'] = d['qc0'] + d['qcf'] + d['m'] + d['cc']

			for k in measures:
				s[k][R] = d[k]

			s['R'] = R

			#print(d)

		for k in measures:
			m[k] = np.mean(s[k])
			m2[k] = np.log2(1+np.mean(s[k]))/(2**N)
			v[k] = np.var(s[k])

		for h in measures:
			row.append(m[h])

		#print(m)
		#print(m2['all'])
	
		table.append(row)

	float_fmt = ".3e"
	table_fmt = 'simple'
	print(tabulate.tabulate(table, headers, floatfmt=float_fmt,
		tablefmt=table_fmt))

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

profile_cpu()
