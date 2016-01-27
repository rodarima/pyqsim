#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys
from qc import *
import tabulate
from T import calc_R

FILE_COM = 'com_qc'
PATH_TEX = '../doc/%s.tex'
PATH_CSV = '../doc/csv/%s.csv'

tabulate.LATEX_ESCAPE_RULES = {}

class QComplexity:
	def __init__(self):
		table_sizes = []
		self.complexity()

	def complexity(self):
		Nmin = 2
		Nmax = 10

		headers = ['n', 'N', 'r', 'R-mean', 'R-var', 'RT-mean', 'RT-var']
		table = []
		for n in range(Nmin, Nmax+1):
			N=2*n
			qs = QSimon({"bits":n})
			config = qs.random_config(n)
			st = qs.run(config)
			qm = QMeasure(st, np.arange(n, n*2))

			cp = CPostprocess(n, qm)
			Rm,Rv,r = cp.profile(10000, timeout=5*60, title='n = %d'%n)
			RTm, RTv = calc_R(n)
			row = [n, N, r, Rm, Rv, RTm, RTv]
			table.append(row)
		print(tabulate.tabulate(table, headers))
		self.save_csv(table, headers, FILE_COM)

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

QComplexity()
