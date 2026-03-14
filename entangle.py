from EMS import *



def rec(Xs,Ys):
	A=min(Xs)
	B=min(Ys)
	C=max(Xs)
	D=max(Ys)
	return(([A,B],[C,B],[C,D],[A,D]))

def sn(M,i,j):
	s=1
	for k in range(i+1,j):
		if M[i][0][0]>= M[k][0][0]:
			s=-s
	return(s)

def sgn(x):
	if x>=0:
		return(1)
	else:
		return(-1)


def entanglements(M):
	E=edges(M,False)
	Ed=[x for x in Edges if abs(x[2])==abs(M[x[0][0]][0][0]-M[x[0][1]][0][0])+abs(M[x[0][0]][0][1]-M[x[0][1]][0][1])]
	return(Ed)



'''


def links(M):
	Ed=entanglements(M)
	links=[]
	for (ij,typ1,mu1) in Ed:
		for (kl,typ2,mu2) in Ed:
			(i,j)=ij
			(k,l)=kl
			if k==j:
				if typ1=="<" and typ2=="<":
						if sgn(mu1)==-sgn(mu2 * sn(M,j,l)):
							links.append((i,j,l),-1)				# sign -1 means that in the spanned square we fill in the lower one (left or right)
						if sgn(mu1)==sgn(mu2 * sn(M,j,l)):
							links.append((i,j,l),1)					# sign 1 means we fill the upper one
					if typ1=="c" and typ2=="c":
						if sgn(mu1)==-sgn(mu2 * sn(M,j,l)):
							links.append((i,j,l),1)	
						if sgn(mu1)==sgn(mu2 * sn(M,j,l)):
							links.append((i,j,l),-1)

def chains(M):
	ents=entanglements(M)
	links=links(M)
	for (ij,typ,m) in ents:
		if ij[0]==0:





def wedge(M,chain):				# input M and a list chain of indices that form a chain
	[A,B]=M[chain[0]][0]
	[C,D]=M[chain[1]][0]
	[E,F]=M[chain[-2]][0]
	[H,I]=M[chain[-1]][0]
	if A<=C and B<=D:
		for (ijl,s) in links(M):
			if ijl[0]==chain[0] and ijl[1]==chain[1] and ijl[2]==chain[3]:
				sign=s
		if sign==-1:
			coch=[i in range(0,len(M)) if E<= M[i][0][0]  and M[i][0][0] <=H and B<= M[i][0][0]  and M[i][0][0] <=D]
			 


'''



def area(M):
	Edges=edges(M,False)
	Ed=[x for x in Edges if abs(x[2])==abs(M[x[0][0]][0][0]-M[x[0][1]][0][0])+abs(M[x[0][0]][0][1]-M[x[0][1]][0][1])]
	Ms=minimals(M)
	Mi=minimali(M)
	Area=[]
	for i in Ms:
		(AB,m)=M[i]
		A=AB[0]
		B=AB[1]
		if abs(m)==A-B+1:
			Area.append(([B,B],[A,B],[A,A],[B,A]))
	for i in Mi:
		(AB,m)=M[i]
		A=AB[0]
		B=AB[1]
		if abs(m)==A-B+1:
			Area.append(([-B,-A],[A,-A],[A,B],[-B,B]))
	for (ij,typ1,mu1) in Ed:
		for (kl,typ2,mu2) in Ed:
			if (ij,typ1,mu1)!=(kl,typ2,mu2):
				(i,j)=ij
				(k,l)=kl
				if k==j:
					[A,B]=M[i][0]
					[C,D]=M[j][0]
					[E,F]=M[l][0]
					if typ1!=typ2:
						Area.append(rec((A,C),(B,D)))
						Area.append(rec((E,C),(F,D)))
					else:
						if typ1=="<":
							if sgn(mu1)==-sgn(mu2 * sn(M,j,l)):
								Area.append(rec((A,C),(B,D)))
								Area.append(rec((E,C),(F,D)))
								Area.append(rec((C,E),(B,D)))
							if sgn(mu1)==sgn(mu2 * sn(M,j,l)):
								#print((A,B,C,D,E,F))
								Area.append(rec((A,C),(B,D)))
								Area.append(rec((E,C),(F,D)))
								Area.append(rec((A,C),(D,F)))
						if typ1=="c":
							if sgn(mu1)==-sgn(mu2 * sn(M,j,l)):
								Area.append(rec((A,C),(B,D)))
								Area.append(rec((E,C),(F,D)))
								Area.append(rec((A,C),(D,F)))	
							if sgn(mu1)==sgn(mu2 * sn(M,j,l)):
								Area.append(rec((A,C),(B,D)))
								Area.append(rec((E,C),(F,D)))
								Area.append(rec((C,E),(B,D)))
	return(Area)







from pa import *

def irarea(E):
	M=conv1(E)
	Amax=max([x[0][0] for x in M])
	Bmax=max([abs(x[0][1]) for x in M])
	points=[]
	for A in range(0,int(Amax)+1):
		for B in range(0,int(Bmax)+1):
			a=A+B+1
			b=A-B+1
			if Is_irred(0,(a,b),E):
				points.append([A,B])
	return(points)


