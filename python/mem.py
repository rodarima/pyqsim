#!/usr/bin/env python3

import numpy as np
import qutip
import scipy.sparse.coo as coo
import time
import sys
import qc
import tabulate

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
		Nmax = 6

		table_qc0 = []
		table_qcf = []
		headers_qc0 = ['N', 'log2_H', 'log2_phi1', 'log2_all_qc0', 'log2_approx_qc0']
		headers_qcf = ['N', 'log2_U', 'log2_qstate', 'log2_all_qcf' ,'log2_approx_qcf']

		for N in range(Nmin, Nmax+1):
			entry = {}
			row_qc0 = []
			row_qcf = []

			#print("N = {}".format(N))
			qs, size_qc0 = self.mem_qc0(N)
			config = qs.random_config(N)
			st, size_qcf = self.mem_qcf(qs, config)
			# Mostrar H y phi1 en qc0
			# U y qstate en qcf
			entry['N']				= N
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

		float_fmt = ".2f"
		table_fmt = 'latex_booktabs'
		print(tabulate.tabulate(table_qc0,
			headers_qc0, floatfmt=float_fmt, tablefmt=table_fmt))
		print(tabulate.tabulate(table_qcf,
			headers_qcf, floatfmt=float_fmt, tablefmt=table_fmt))

QMemProfiler()
