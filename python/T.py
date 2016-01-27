def I(i,n):
	return 1 - 2**(i-n+1)

def _I(i,n):
	return 2**(i-n+1)

def gb(m, r, q):
	if r>m: return 0
	num = 1
	for i in range(m-r+1, m+1):
		num*=(1-q**i)
	denom=1
	for i in range(1, r+1):
		denom*=(1-q**i)
	return num//denom

def T(i, p, n):
	s = 1.0
	for j in range(0, n-1):
		s *= I(j,n)

	d = (p-n+1)
	if(d==0): return s

	s *= 2**((-n+1)*d)
	#l = [0] * d
	#li = [0] * (d*(n-2)+1)
	#sd = 0.0
	#end = False
	#while True:
	#	#print(l)

	#	pd = 1.0
	#	#for j in range(d):
	#	#	_i = 2**l[j]
	#	#	#_i = _I(l[j], n)
	#	#	#print(_i)
	#	#	pd *= _i
	#	#sd += pd

	#	si = 0
	#	for j in range(d):
	#		si += l[j]
	#	#print(si)
	#	li[si]+=1
	#	sd += 2**si

	#	j = d-1
	#	while True:
	#		if l[j] < (n-2):
	#			l[j]+=1
	#			e = l[j]
	#			for k in range(j+1, d):
	#				l[k] = e
	#			break

	#		if j==0:
	#			end = True
	#			break

	#		j = j-1

	#	if end: break

	#s *= sd
	#print(sd, (p-1, p-n+1, 2), gb(p-1,p-n+1,2))#, sum(li), li)

	s *= gb(p-1,p-n+1,2)

	return s

def test():
	for n in range(2, 12):
		E = 0
		E2 = 0
		for r in range(n-1, 100):
			p = T(0, r, n)
			E += r * p
			E2 += r**2 * p
		Var = E2 - E**2
		print("N={}\tE[R]={:.4f}\tE/n={:.3f}\tVar={:.3f}".format(n, E, E/n, Var))

def calc_R(n):
	E = 0
	E2 = 0
	for r in range(n-1, 100):
		p = T(0, r, n)
		E += r * p
		E2 += r**2 * p
	Var = E2 - E**2
	return (E, Var)

