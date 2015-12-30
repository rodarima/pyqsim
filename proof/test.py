def I(i,n):
	return 1-2**(i-n+1)

def _I(i,n):
	return 1-I(i,n)

def T(i,p,n):
	if (p == 0) and (i != n-1): return 0
	elif (i == n-1) and (p != 0): return 0
	elif p == 0: return 1
	else: return _I(i,n)*T(i,p-1,n)+I(i,n)*T(i+1,p-1,n)

def test():
	for n in range(2, 10):
		E = 0
		for r in range(n-1, 30):
			E += r*T(0, r, n)
		print("N={}\tE[R]={}".format(n, E))

test()
