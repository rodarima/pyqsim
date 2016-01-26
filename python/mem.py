#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys
import qc
import tabulate

FILE_QC0 = 'table_qc0'
FILE_QCF = 'table_qcf'
PATH_TEX = '../doc/%s.tex'
PATH_CSV = '../doc/csv/%s.csv'

tabulate.LATEX_ESCAPE_RULES = {}

class QMemProfiler:
	def __init__(self):
		table_sizes = []
		self.profile_mem()

	def mem_qc0(self, N):
		qs = qc.QSimon({"bits":N})
		sizes = qs.sz
		#print('Memory QC0: {}'.format(sizes))
		return qs, sizes

	def mem_qcf(self, qs, config):
		st = qs.run(config)
		sizes = qs.sz
		#print('Memory QCf: {}'.format(sizes))
		return st, sizes

	def profile_mem(self):
		Nmin = 2
		Nmax = 10

		table_qc0 = []
		table_qcf = []
		headers_qc0 = ['n', 'N', 'log2_H', 'log2_phi1', 'log2_all_qc0', 'log2_approx_qc0']
		headers_qcf = ['n', 'N', 'log2_H', 'log2_U', 'log2_qstate', 'log2_all_qcf' ,'log2_approx_qcf']
		tex_headers_qc0 = [
			'$n$',
			'$N$',
			'$\log_2 S(H)$',
			'$\log_2 S(\ket{\phi_1})$',
			'$\log_2 S_T$',
			'$\log_2 S_T\'$'
		]
		tex_headers_qcf = [
			'$n$',
			'$N$',
			'$\log_2 S(H)$',
			'$\log_2 S(U)$',
			'$\log_2 S(\ket{\phi_3})$',
			'$\log_2 S_T$'
			,'$\log_2 S_T\'$'
		]

		for N in range(Nmin, Nmax+1):
			entry = {}
			row_qc0 = []
			row_qcf = []

			print("N = {}".format(N))
			qs, size_qc0 = self.mem_qc0(N)
			config = qs.random_config(N)
			st, size_qcf = self.mem_qcf(qs, config)
			# Mostrar H y phi1 en qc0
			# U y qstate en qcf
			entry['N']				= 2*N
			entry['n']				= N
			entry['H']				= size_qc0['H']
			entry['phi1']			= size_qc0['phi1']
			entry['all_qc0']		= size_qc0['all']
			entry['log2_H']			= np.log2(size_qc0['H'])
			entry['log2_phi1']		= np.log2(size_qc0['phi1'])
			entry['log2_all_qc0']	= np.log2(size_qc0['all'])
			entry['log2_approx_qc0'] = np.log2(2**(2*N+3) + 2**(2*N+2))
			# El estado phi1 puede ser machacado?
			entry['U']				= size_qcf['U']
			entry['qstate']			= size_qcf['qstate']
			entry['all_qcf']		= size_qcf['U'] + size_qcf['H'] + size_qcf['qstate']
			entry['log2_U']			= np.log2(size_qcf['U'])
			entry['log2_qstate']	= np.log2(size_qcf['qstate'])
			entry['log2_all_qcf'] = np.log2(
				size_qcf['U'] + size_qcf['H'] + size_qcf['qstate'])
			entry['log2_approx_qcf'] = np.log2(
#				2**(2*N+3) + 2**(2*N+2) +	#U
#				2**(2*N+3) +				#H
#				2**(2*N+2) + 2**(2*N+1)		#qstate
				2**(2*N+3) + 2**(2*N+2) +
				2**(2*N+3) +
				2**(2*N+2) + 2**(2*N+1)
			)

			for h in headers_qc0:
				row_qc0.append(entry[h])
			for h in headers_qcf:
				row_qcf.append(entry[h])

			table_qc0.append(row_qc0)
			table_qcf.append(row_qcf)

		self.save_tex(table_qc0, tex_headers_qc0, FILE_QC0)
		self.save_tex(table_qcf, tex_headers_qcf, FILE_QCF)

		self.save_csv(table_qc0, headers_qc0, FILE_QC0)
		self.save_csv(table_qcf, headers_qcf, FILE_QCF)

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

QMemProfiler()
