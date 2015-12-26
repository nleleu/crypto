#!/usr/bin/python
import math
import sys
import decimal
import numbers
def is_prime(n):
    if n % 2 == 0 and n > 2: 
        return False
    i = 3
    g = decimal.Decimal(str(n))
    while i <  int(g.sqrt()) + 1:
        if n % i == 0:
            return False
	i += 2
    return True

def crack(n):
	i = 1
	g = decimal.Decimal(str(n))
	while i < int(g.sqrt()):
		if is_prime(i):
			di = decimal.Decimal(str(i))
			result = g/di
			if result._isinteger() and is_prime(n/i):
				print "test de "+str(i)+" et "+str(n/i)
				p = i
				q = n/i
				res = bezout(n_global, (p - 1)*(q - 1))
				u = res[0]
				v = res[1]
				print "p = "+str(p)
				print "q = "+str(q)
				print "u = "+str(u)
				print "v = "+str(v)
				return
		i +=1				
		

def bezout(a, b):
	r = a
	rp = b
	u = 1
	v = 0
	up = 0
	vp = 1
	while rp != 0:
		q = int(r/rp)
		rs = r
		us = u
		vs = v
		r = rp
		u = up
		v = vp
		rp = rs -q * rp
		up = us -q*up
		vp = vs - q * vp			
	return (u,v)

n_global = int(sys.argv[1])
print n_global
crack(n_global)

