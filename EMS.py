from pa import *
import math


def si(x):						#sign function
	if x>=0:
		return(1)
	else:
		return(-1)



def sign_con(M):		#sign-condition for EMS
	#print('hi')
	s=0
	for i in range(len(M)):
		b=M[i][0][0]-M[i][0][1]+1
		mu=M[i][1]
	#	print(mu)
		s=s+math.floor(mu/2)
		for j in range(i):
			bj=M[j][1]
			s=(s+mu*(bj-1))
	#	print(s)
	if s%2==0:	
		return(True)
	else:
		return(False)


def check(M):					#checks if M is actually an extended multisegment
	for (AB,m) in M:
		if (AB[0]-AB[1]+1-m)%2!=0:
			return(False)
		if m>abs(AB[0]-AB[1]+1):
			return(False)
	for i in range(0,len(M)-1):
		[A,B]=M[i][0]
		for j in range(i+1,len(M)):
			[C,D]=M[j][0]
			if C<A and D<B:
				return(False)
	#if not sign_con(M):
	#	return(False)
	else:
		return(True)




def conv1(M):						#convert extended multisegment M to other parametrisation
	E=[]
	s=1
	for ([A,B],l,e) in M:
		b=A-B+1
		m=s*e*(b-2*l)
		s=s*(-1)**(b-1)
		E.append(([A,B],m))
	return(E)	



def conv2(E):						#convert it back to Atobe parametrisation
	E=sort(E)
	M = []
	s = 1
	for ([A, B], m) in E:
		b = A - B + 1
		m_prime = m / s 
		e_sign = 1 if m_prime >= 0 else -1
		m_prime_abs = abs(m_prime)
		l_numerator = b - m_prime_abs
		l = int(l_numerator / 2)
		e = e_sign 
		M.append(([A, B], l, e))
		s = s * ((-1) ** (b - 1))
	return M	





def swap(M,i,check=False):					#exchanges the order of two consecutive elements of extended multisegment M where the first one (left) has index i
	(AB,m)=M[i]
	(CD,n)=M[i+1]
	[A,B]=AB 
	[C,D]=CD
	j=len(M)
	if A<C and B<D:
		if check:
			return(M,False)
		else:
			return(M)
	if A>=C and B<=D:
		E=[]
		for r in range(0,j):
			if r==i:
				E.append(([C,D],n))
				mu=2*n-m
				E.append(([A,B],mu))
			if r==i+1:
				e=1
			else:
				if r!=i and r!=i+1:
					E.append(M[r])
		if check:
			return(E,True)
		else:
			return(E)
	if A<=C and B>=D:
		E=[]
		for r in range(0,j):
			if r==i:
				mu=2*m-n
				E.append(([C,D],mu))
				E.append(([A,B],m))
			if r==i+1:
				e=1
			else:
				if r!=i and r!=i+1:
					E.append(M[r])
		if check:
			return(E,True)
		else:
			return(E)




def betas(M):
	betas=[]
	j=0
	while j<len(M)-1:
		[A,B]=M[j][0]
		[C,D]=M[j+1][0]
		if (D<B and A<=C) or (D==B and A<C):
			betas.append(j)
			return(betas)
		j=j+1
	return(betas)

def sort(A):				#brings extended multisegment into standard form
	E=A[:]
	bet=betas(E)
	while len(bet)!=0:
		x=bet[0]
		E=swap(E,x)
		bet=betas(E)
	return(E)

def ui(M,i):					#union intersection or the inverse at index i 
	(AB,m)=M[i]
	(CD,n)=M[i+1]
	[A,B]=AB 
	[C,D]=CD
	j=len(M)
	if abs(m-n)<abs(A-C)+abs(B-D):
		return(M)
	elif B<=D:
		E=[]
		for r in range(0,j):
			if r==i:
				E.append(([C,B],m-si(n-m)*(C-A)))
				E.append(([A,D],n-si(n-m)*(C-A)))
			elif r==i+1:
				e=1
			else:
				E.append(M[r])
		if check(E):
			return(E)
		else:
			return(M)
	elif A<=C and B>=D:
		M=swap(M,i)
		return(ui(M,i))


def breakup(M,i,x):						#the inverse of union intersection applied to a single extended segment
	mu=M[i][1]
	[A,B]=M[i][0]
	if abs(mu)!= A-B+1:
		return(M)
	else:
		if not (A >x and B<=x):
			return(M)
		else:
			E=[]
			for j in range(0,i):
				E.append(M[j])
			E.append(([x,B],mu-si(mu)*(A-x)))
			E.append(([A,x+1],-si(mu)*(A-x)))
			for j in range(i+1,len(M)):
				if (int(2*A))%2==0:
					E.append((M[j][0],-M[j][1]))
				else:
					E.append((M[j][0],M[j][1]))
	if check(E):
		return(E)
	else:
		return(M)




def diff(M,i,j):			#M must be standard, computes the difference between extended segments i and j 
	if j==i+1:
		return(M[j][1]-M[i][1])
	else:
		A,B=M[i][0]
		m=M[j][1]-M[j-1][1]
		C,D=M[j-1][0]
		if A>=C and B<=D:
			return(m-diff(M,i,j-1))
		else:
			return(m+diff(M,i,j-1))


def strict(M,i,j):  #M must be in standard form for this function, checks wheter i and j are connected
	[A,B]=M[i][0]
	[C,D]=M[j][0]
	if A<C and B<D:
		for k in range(i+1,j):
			if A<M[k][0][0] and M[k][0][0]<C and B<M[k][0][1] and M[k][0][1]<D:
				return(False) 
	if A>=C and B<=D:
		for k in range(i+1,j):
			if A>=M[k][0][0] and M[k][0][0]>=C and B<=M[k][0][1] and M[k][0][1]<=D:
				return(False)  
	return(True)



def join(M,i,j):			#M must be standard to call this function, which makes segment i and j adjacent
	if j==i+1:
		return(M,i)
	elif M[i][0][0]>=M[i+1][0][0] and M[i][0][1]<=M[i+1][0][1]:
		return(join(swap(M,i),i+1,j))
	elif M[j-1][0][0]>=M[j][0][0] and M[j-1][0][1]<=M[j][0][1]:
		return(join(swap(M,j-1),i,j-1))


############################################################################





def edges(M, rest=False):   #returns the edge pairs where M is an extended multi-segment
	edges=[]
	M=sort(M)
	for i in range(0,len(M)):
		for j in range(i+1,len(M)):
			if strict(M,i,j):
				if M[i][0][0]>=M[j][0][0]:
					typ="c"
				else:
					typ="<"
				if not rest:
					edges.append(((i,j),typ,diff(M,i,j)))
				else:
					edges.append(((i,j),typ,abs(M[i][0][0]-M[j][0][0])+ abs(M[i][0][1]-M[j][0][1]) - abs(diff(M,i,j))))
	return(edges)


def minimals(M):
	m_values=[]
	for i in range(0,len(M)):
		minimal=True
		for j in range(i+1,len(M)):
			if M[i][0][0]>=M[j][0][0] and M[i][0][1]<=M[j][0][1]:
				minimal=False
		if minimal:
			m_values.append(i)
	return(m_values)


def minimali(M):
	m_values1=[]
	for i in range(0,len(M)):
		minimal1=True
		for j in range(0,i):
			if M[j][0][0]<M[i][0][0] and M[j][0][1]<M[i][0][1]:
				minimal1=False
		if minimal1:
			m_values1.append(i)
	return(m_values1)


def frp(A,B,M):
	E=conv2(M)
	a=int(A+B+1)
	b=int(A-B+1)
	return(FRP((a,b),E))




###############################################################################
###############################################################################


from itertools import product

def swap_forward(M, i):	#makes the segment at index i the last amongst the ones with the same B-value
	M=sort(M)
	E = list(M)
	pos = i
	while pos < len(E) - 1:
		(AB, m) = E[pos + 1]
		A, B = AB
		(CD, n) = E[pos]
		C, D = CD
		if D > B:
			break
		E, changed = swap(E, pos, check=True)
		if not changed:
			break
		pos += 1
	return (E,pos)

def maxi(M):
	M=sort(M)
	ed=edges(M,True)
	holds=True
	for i in range(len(M)):
		(AB,m)=M[i]
		#print(abs(m),AB[0]-AB[1]+1)
		if int(abs(m))==int(AB[0]-AB[1]+1) and abs(m)!=1:
			M,pos=swap_forward(M,i)
			E=breakup(M,pos,M[pos][0][1])
			if E!=M:
				return(maxi(sort(E)))
	for E in ed:
		if E[1]=="c" and E[2]==0:
			AB=M[E[0][0]][0]
			CD=M[E[0][1]][0]
			if AB[0]!=CD[0] and AB[1]!=CD[1]: 
				M,i=join(M,E[0][0],E[0][1])
				E=ui(M,i)
				if E!=M:
					return(maxi(sort(E)))
	return(sort(M))


