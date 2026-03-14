# To load a file "packet.sage", use:
# load("packet.sage")


def Notation():
    print("Notation:")
    print("m = ([x_1,y_1], [x_2,y_2], ... ,[x_r,y_r]) is a multi-segment.")
    print("For example, m = ([1/2,-1/2],) is the Steinberg of GL(2).") 
    print("")
    print("T = ([x_1,e_1], ... ,[x_t,e_t]) is a tempered representation of G")
    print(r"with L-parameter \phi = S_{2x_1+1} + ... + S_{2x_t+1}")
    print(r"and character \ep(S_{2x_i+1}) = e_i.")
    print("For example, T = ([1/2,-1],[3/2,+1],[5/2,-1]) is a supercuspidal representation of SO(13).")
    print("")
    print("If T is trivial, write T = ([1,1],).")
    print("")
    print("psi = ((a_1,b_1), ... ,(a_r,b_r)) is an A-parameter,")
    print(r"which means that S_{a_1} \otimes S_{b_1} + ... + S_{a_r} \otimes S_{b_r}.")
    print("")
    print("E = (([A_1,B_1],l_1,eta_1), ... ,([A_m,B_m],l_m,eta_m)) is an extended multi-segment.")

#psi = ((2,2),(5,3))
#psi1 = ((1,6),(1,2),(4,1))
psi2 = ((51,31),(31,45),(13,5))
#psi3 = ((1, 5), (3, 3), (4, 4))
#psi4 = ((2, 5), (3, 2), (5, 4))
#E = (([0,0],0,-1),([2,0],0,-1),([2,2],0,-1))
#(a,b) = (4,2)

phi = ((1,1),(3,1),(5,1))
phi1 = ((1, 1), (3, 1), (5, 1), (7, 1), (9, 1))
phi2 = ((2, 1), (4, 1), (6, 1), (8, 1))
m = ([0,-3],[2,-3])
m1 = ([-2,-3],[-1,-1])
m2 = ([-1/2,-5/2],)

E0 = (([5/2, -3/2], 2, -1), ([3/2, 1/2], 0, -1), ([7/2, 1/2], 2, 1))
psi = ((2, 5), (3, 2), (5, 4))
E4 = (([5,0],1,-1), ([7,1],3,+1))


###########################################################
# import. 
###########################################################

import numpy as np

###########################################################
# The best matching function for the left derivatives. 
###########################################################

def LBM(x,m): 
    """ 
    Left best matching functions.       
    Input: (x,m)
        x (real number)
        m (multi-segment)
    Output: (A,B,A0,B0,A1,B1)
        A  ( = A_{x-1} )
        B  ( = A_{x} )
        A0 ( = A_{x-1}**0 )
        B0 ( = A_{x}**0 )
        A1 ( = A_{x-1}**c )
        B1 ( = A_{x}**c )
    """
    A = tuple(filter(lambda z: z[0] == x-1, list(m)))
    B = tuple(filter(lambda z: z[0] == x, list(m)))
    A0 = []; B0 = []; C0 = list(B)
    for i in range(len(A)): 
        for j in range(len(C0)): 
            if A[len(A)-i-1][1] < C0[j][1]:
                A0.append(A[len(A)-i-1])
                B0.append(C0[j])
                del C0[j]
                break
    A0 = tuple(A0); B0 = tuple(B0)
    A1 = list(A); B1 = list(B); A2 = list(A0); B2 = list(B0)
    for j in range(len(A2)): 
        for i in range(len(A1)):
            if A1[i] == A2[j]:
                del A1[i]
                break
    for j in range(len(B2)): 
        for i in range(len(B1)):
            if B1[i] == B2[j]:  
                del B1[i]
                break
    A1 = tuple(A1); B1 = tuple(B1)
    return A, B, A0, B0, A1, B1


##################################################################
# The left derivatives and socles for GL(n). 
##################################################################

def LD(x,m):
    """
    Left highest derivative.
    Input: (x,m)
        x (real number)
        m (multi-segment)
    Output: (k,n)
        k (non-negative integer)
        n (multi-segment)
            L_{x}**{max}(m) = L_{x}**{(k)}(m) = n
    """
    A, B, A0, B0, A1, B1 = LBM(x,m)
    C = list(m)
    for j in range(len(B1)):
        for i in range(len(C)): 
            if C[i] == B1[j]: 
                C[i] = [x-1,C[i][1]]
                break
    n = tuple(sorted(filter( lambda z: z[0] >= z[1], C)))
    return len(B1), n

def LS(x,m):
    """
    Left socle.
    Input: (x,m)
        x (real number)
        m (multi-segment)
    Output: n
        n (multi-segment)
            soc( ||**{x} \times m ) = n
    """
    A, B, A0, B0, A1, B1 = LBM(x,m)
    C = list(m)
    if A1 == ():
        C.append([x,x])
    else: 
        for i in range(len(C)):
            if C[i] == A1[0]:
                C[i] = [x,C[i][1]]
                break
    n = tuple(sorted(C))
    return n 


###########################################################
# The best matching function for the right derivatives. 
###########################################################

def RBM(y,m): 
    """ 
    Right best matching functions.       
    Input: (y,m)
        y (real number)
        m (multi-segment)
    Output: (A,B,A0,B0,A1,B1)
        A  ( = B_{y+1} )
        B  ( = B_{y} )
        A0 ( = B_{y+1}**0 )
        B0 ( = B_{y}**0 )
        A1 ( = B_{y+1}**c )
        B1 ( = B_{y}**c )
    """
    A = tuple(filter(lambda z: z[1] == y+1, list(m)))
    B = tuple(filter(lambda z: z[1] == y, list(m)))
    A0 = []; B0 = []; C0 = list(B)
    for i in range(len(A)): 
        for j in range(len(C0)): 
            if A[i][0] > C0[len(C0)-j-1][0]:
                A0.append(A[i])
                B0.append(C0[len(C0)-j-1])
                del C0[len(C0)-j-1]
                break
    A0 = tuple(A0); B0 = tuple(B0)
    A1 = list(A); B1 = list(B); A2 = list(A0); B2 = list(B0)
    for j in range(len(A2)): 
        for i in range(len(A1)):
            if A1[i] == A2[j]:
                del A1[i]
                break
    for j in range(len(B2)): 
        for i in range(len(B1)):
            if B1[i] == B2[j]:  
                del B1[i]
                break
    A1 = tuple(A1); B1 = tuple(B1)
    return A, B, A0, B0, A1, B1


###########################################################
# The right derivatives and socles for GL(n). 
###########################################################

def RD(y,m):
    """
    Right highest derivative.
    Input: (y,m)
        y (real number)
        m (multi-segment)
    Output: (k,n)
        k (non-negative integer)
        n (multi-segment)
            R_{x}**{max}(m) = R_{x}**{(k)}(m) = n
    """
    A, B, A0, B0, A1, B1 = RBM(y,m)
    C = list(m)
    for j in range(len(B1)):
        for i in range(len(C)): 
            if C[i] == B1[j]: 
                C[i] = [C[i][0],y+1]
                break
    n = tuple(sorted(filter(lambda z: z[0] >= z[1], C)))
    return len(B1), n

def RS(y,m):
    """
    Left socle.
    Input: (x,m)
        x (real number)
        m (multi-segment)
    Output: n
        n (multi-segment)
            soc( m \times ||**{x} ) = n
    """
    A, B, A0, B0, A1, B1 = RBM(y,m)
    C = list(m)
    if A1 == ():
        C.append([y,y])
    else: 
        for i in range(len(C)):
            if C[i] == A1[len(A1)-1]:
                C[i] = [C[i][0],y]
                break
    n = tuple(sorted(C))
    return n


###########################################################
# The x-derivatives for G (x != 0).
###########################################################

def D(x,m,T):
    """
    Highest x-derivative.
    Input: (x,m,T)
        x (nonzero real number)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (k,dm,dT)
        k (non-negative integer)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            D_{x}**{max}(m,T) = D_{x}**{(k)}(m,T) = (dm,dT)
    """
    if x == 0:
        print("We cannot compute the 0-derivatives.")
    if x < 0: 
        k,n = LD(x,m); DT = T
        return k,n,DT
    if x > 0:
        T0 = list(T)
        if 2*x % 2 == 1:
            T0.append([-1/2,+1]) 
        t = m.count([x-1,-x])
        n = list(filter(lambda z: z != [x-1,-x], m))
        e = (-1)**t  
        A, B, A0, B0, A1, B1 = LBM(x,n)
        C, D, C0, D0, C1, D1 = RBM(-x,n)
        for j in range(len(B1)):        
            for i in range(len(n)):     
                if n[i] == B1[j]:       
                    n[i] = [x-1,n[i][1]]
                    break               
        if T0.count([x,1])+T0.count([x,-1]) == 0:
            m1 = T0.count([x,1])+T0.count([x,-1]) 
            m0 = T0.count([x-1,1])+T0.count([x-1,-1])                           
            dt = t            
        if T0.count([x,1])*T0.count([x-1,-e]) > 0:
            m1 = T0.count([x,1]) -1                
            m0 = T0.count([x-1,-e]) -1
            for i in range(max(m1-len(A1),0)): 
                T0.remove([x,1])           
                T0.append([x-1,-e])
            if max(m1-len(A1),0) % 2 == 0:
                dt = t
            else:                             
                T0.remove([x,1])
                T0.remove([x-1,-e])
                dt = t+1 
        if T0.count([x,-1])*T0.count([x-1,e]) > 0:
            m1 = T0.count([x,-1]) -1                
            m0 = T0.count([x-1,e]) -1
            for i in range(max(m1-len(A1),0)): 
                T0.remove([x,-1])           
                T0.append([x-1,e])
            if max(m1-len(A1),0) % 2 == 0:
                dt = t
            else:                             
                T0.remove([x,-1])
                T0.remove([x-1,e])
                dt = t+1 
        if T0.count([x,1])*T0.count([x-1,e]) > 0:   
            m1 = T0.count([x,1])                     
            m0 = T0.count([x-1,e])                     
            for i in range(max(m1-len(A1),0)):           
                T0.remove([x,1])                     
                T0.append([x-1,e])                     
            if max(m1-len(A1),0) % 2 == 0 or t == 0:     
                dt = t           
            else:                                      
                T0.append([x,1])                                     
                T0.append([x-1,e])                  
                dt = t-1      
        if T0.count([x,-1])*T0.count([x-1,-e]) > 0: 
            m1 = T0.count([x,-1])                    
            m0 = T0.count([x-1,-e])                    
            for i in range(max(m1-len(A1),0)):          
                T0.remove([x,-1])                   
                T0.append([x-1,-e])                   
            if max(m1-len(A1),0) % 2 == 0 or t == 0:    
                dt = t         
            else:                                       
                T0.append([x,-1])                   
                T0.append([x-1,-e])                   
                dt = t-1        
        if T0.count([x-1,1])+T0.count([x-1,-1]) == 0 and T0.count([x,1]) > 0:
            m1 = T0.count([x,1]); m0 = 0                                       
            for i in range(max(m1-len(A1),0)):                                    
                T0.remove([x,1])                                              
                T0.append([x-1,e])                                              
            if max(m1-len(A1),0) % 2 == 0 or t == 0:                              
                dt = t                                    
            else:                                                                 
                T0.append([x,1])                                              
                T0.append([x-1,e])                                              
                dt = t-1                                  
        if T0.count([x-1,1])+T0.count([x-1,-1]) == 0 and T0.count([x,-1]) > 0:
            m1 = T0.count([x,-1]); m0 = 0                                       
            for i in range(max(m1-len(A1),0)):                                     
                T0.remove([x,-1])                                              
                T0.append([x-1,-e])                                              
            if max(m1-len(A1),0) % 2 == 0 or t == 0:                               
                dt = t                                    
            else:                                                                  
                T0.append([x,-1])                                              
                T0.append([x-1,-e])                                              
                dt = t-1      
        for j in range(len(D1)-m0-max(len(A1)-m1,0)):                              
            for i in range(len(n)):                                                
                if n[i] == D1[j]:                                                  
                    n[i] = [n[i][0],-(x-1)]                                        
                    break                                                          
        for j in range(dt):                                                        
            n.append([x-1,-x]) 
        dm = tuple(sorted(filter(lambda z: z[0] >= z[1], n)))
        dk = len(B1)+max(m1+max(len(D1)-m0,0)-len(A1),0)
        DT = tuple(sorted(filter(lambda z: z[0] != -1/2, T0)))
        return dk, dm, DT


###########################################################
# The x-derivatives: The case where k <= 1.
###########################################################

def D1(x,m,T):
    """
    Highest x-derivative: The case where k <= 1
    Input: (x,m,T)
        x (nonzero real number)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (dm,dT)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            D_{x}**{(1)}(m,T) = (dm,dT)
    """
    T0 = list(T)
    t = m.count([x-1,-x])
    n = list(filter(lambda z: z != [x-1,-x], m))
    A, B, A0, B0, A1, B1 = LBM(x,n)
    if len(B1) > 0: 
        for i in range(len(n)):
            if n[i] == B1[0]:
                n[i] = [x-1,n[i][1]]
                break 
    else:
        e = (-1)**t
        if 2*x % 2 == 1:
            T0.append([-1/2,+1])
        if T0.count([x,1])+T0.count([x,-1]) == 0:
            c = 0
            m1 = 0
            m0 = T0.count([x-1,1])+T0.count([x-1,-1])                                       
        elif T0.count([x,1])*T0.count([x-1,-e]) > 0:
            c = 1
            m1 = T0.count([x,1]) -1                
            m0 = T0.count([x-1,-e]) -1
        elif T0.count([x,-1])*T0.count([x-1,e]) > 0:
            c = 2
            m1 = T0.count([x,-1]) -1                
            m0 = T0.count([x-1,e]) -1
        elif T0.count([x,1]) > 0 and T0.count([x-1,-e]) == 0:   
            c = 3
            m1 = T0.count([x,1])                     
            m0 = T0.count([x-1,e])                     
        elif T0.count([x,-1]) > 0 and T0.count([x-1,e]) == 0: 
            c = 4
            m1 = T0.count([x,-1])                    
            m0 = T0.count([x-1,-e])                    
        if m1 > len(A1): 
            if c == 1:
                T0.remove([x,1])
                T0.remove([x,1])
                t = t+1
            elif c == 2:
                T0.remove([x,-1])
                T0.remove([x,-1])
                t = t+1
            elif t > 0 and c == 3:
                T0.append([x-1,e])
                T0.append([x-1,e])
                t = t-1
            elif t > 0 and c == 4:
                T0.append([x-1,-e])
                T0.append([x-1,-e])
                t = t-1
            elif t == 0 and c == 3:
                T0.remove([x,1])           
                T0.append([x-1,1])
                t = 0
            elif t == 0 and c == 4:
                T0.remove([x,-1])           
                T0.append([x-1,-1])
                t = 0
        else: 
            C, D, C0, D0, C1, D1 = RBM(-x,n)
            if m1+len(D1)-m0-len(A1) > 0:                              
                for i in range(len(n)):                                                
                    if n[i] == D1[0]:                                                  
                        n[i] = [n[i][0],-(x-1)]                                        
                        break                                                                      
            else: 
                n = []; T0 = []; t = 0
    for j in range(t):                                                        
        n.append([x-1,-x])                                                     
    dm = tuple(sorted(filter(lambda z: z[0] >= z[1], n)))
    DT = tuple(sorted(filter(lambda z: z[0] != -1/2, T0)))
    return dm, DT


###########################################################
# The x-socles for G (x != 0).
###########################################################

def S(x,m,T):              
    """
    x-socle.
    Input: (x,m,T)
        x (nonzero real number)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (dm,dT)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            soc( ||**{x} \rtimes (m,T) ) = (dm,dT)
    """
    t = m.count([x-1,-x])
    n = list(filter(lambda z: z != [x-1,-x], m))
    T0 = list(T); e = (-1)**t
    if 2*x % 2 == 1:
        T0.append([-1/2,+1])                                                             
    if T0.count([x,1])*T0.count([x-1,-e]) + T0.count([x,-1])*T0.count([x-1,e]) > 0:
        m1 = T0.count([x,1])+T0.count([x,-1]) -1
        m0 = T0.count([x-1,1])+T0.count([x-1,-1]) -1
    else:                                           
        m1 = T0.count([x,1])+T0.count([x,-1]) 
        m0 = T0.count([x-1,1])+T0.count([x-1,-1]) 
    A, B, A0, B0, A1, B1 = LBM(x,n) 
    C, D, C0, D0, C1, D1 = RBM(-x,n)
    if x < 0:                      
        if len(A1) > 0:                  
            for i in range(len(n)):      
                if n[i] == A1[0]:
                    n[i] = [x,n[i][1]]
                    break  
        else:              
            n.append([x,x])    
        dt = t                                                    
    if x > 0:                                                                              
        if m1 + max(len(D1)-m0,0) < len(A1):                                               
            for i in range(len(n)):                                                        
                if n[i] == A1[0]:             
                    n[i] = [x,n[i][1]]                
                    break                           
            dt = t          
        elif len(D1) < m0 and m1 >= len(A1):         
            if T0.count([x-1,1]) > 0:       
                if T0.count([x-1,1]) == m0 and t == 0:
                    dt = 0 
                    T0.append([x,e])
                    T0.remove([x-1,1])                                
                elif T0.count([x-1,1]) == m0 and t > 0: 
                    dt = t-1                 
                    T0.append([x,e])     
                    T0.append([x,e])     
                else:                        
                    dt = t+1                 
                    T0.remove([x-1,1])     
                    T0.remove([x-1,1])                                
            if T0.count([x-1,-1]) > 0:                                 
                if T0.count([x-1,-1]) == m0 and t == 0:
                    dt = 0                                            
                    T0.remove([x-1,-1])                               
                    T0.append([x,-e])                               
                elif T0.count([x-1,-1]) == m0 and t > 0:                    
                    dt = t-1                                            
                    T0.append([x,-e])                               
                    T0.append([x,-e])                               
                else:                                                   
                    dt = t+1                 
                    T0.remove([x-1,-1])    
                    T0.remove([x-1,-1])    
        elif len(D1) >= m0 and m1+len(D1)-m0 >= len(A1) and len(C1) > 0:
            for i in range(len(n)):                                     
                if n[i] == C1[len(C1)-1]:                                       
                    n[i] = [n[i][0],-x]                                 
                    break                                               
            dt = t                              
        elif len(D1) >= m0 and m1+len(D1)-m0 >= len(A1) and len(C1) == 0:                                                  
            n.append([-x,-x])                                           
            dt = t
    for j in range(dt):
        n.append([x-1,-x])
    dm = tuple(sorted(filter(lambda z: z[0] >= z[1], n)))
    DT = tuple(sorted(filter(lambda z: z[0] != -1/2, T0)))
    return dm, DT


###########################################################
# The [0,-1]-derivatives and socles for G.
###########################################################

def D00(m,T):
    """
    Highest [0,-1]-derivative.
    Input: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
            (m,T) is ||**{-1} reduced.
    Output: (k,dm,dT)
        k (non-negative integer)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            D_{[0,-1]}**{max}(m,T) = D_{[0,-1]}**{(k)}(m,T) = (dm,dT)
    """
    if D(-1,m,T)[0] > 0:
        print("L(m; T) has a nontrivial (-1)-derivative.")
    else:
        k0,m0 = LD(0,m)
        k1,m1 = LD(-1,m0)
        for i in range(k0-k1):
            m1 = LS(0,m1)
        return k1,m1,T

def S00(m,T): 
    """
    [0,-1]-socle.
    Input: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
            (m,T) is ||**{-1} reduced.
    Output: (dm,dT)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            soc( Delta[0,-1] \rtimes (m,T) ) = (dm,dT)
    """
    if D(-1,m,T)[0] > 0:
        print("L(m; T) has a nontrivial (-1)-derivative.")
    else: 
        k0,m0 = LD(0,m)
        m1 = LS(-1,m0)
        for i in range(k0+1):
            m1 = LS(0,m1)
        return m1,T


###########################################################
# The [0,1]-derivatives for G: A special case.
###########################################################

def D01A(s,t,T):                                                
    """
    Highest [0,1]-derivative: A special case. 
    Input: (s,t,T)
        s (non-negative integer)
        t (non-negative integer)
        T (tempered L-parameter)
    Output (k, ds, dt, dT)
        k (non-negative integer)
        ds (non-negative integer)
        dt (non-negative integer)                
        dT (tempered L-parameter)
            D_{[0,1]}**{(k)}( L((||**1)**{s} \times [0,-1]**{t} \times T) )
            = L((||**1)**{ds} \times [0,-1]**{dt} \times dT)
    """
    n = []
    for j in range(s):                                        
        n.append([-1,-1])                                                        
    for j in range(t):                                        
        n.append([0,-1])
    m = tuple(n)                                      
    if D(1,m,T)[0] > 0:                                       
        print("L([-1,-1]**s,[0,-1]**t; T) has a nontrivial 1-derivative.")
    else:    
        e = (-1)**t; T0 = list(T)
        if T.count([1,1]) > 0:
            if (T.count([0,-e])-s) % 2 == 1: 
                if t % 2 == 0: 
                    dk = t; ds = s; dt = 0 
                    DT = tuple(sorted(T0))
                else: 
                    dk = t; ds = s+1; dt = 0
                    T0.append([0,-e])
                    T0.remove([1,1])
                    DT = tuple(sorted(T0)) 
            else: 
                if t % 2 == 1:
                    dk = t+1; ds = s; dt = 0
                    T0.remove([0,-e])
                    T0.remove([1,1])
                    DT = tuple(sorted(T0)) 
                elif s > 0:
                    dk = t+1; ds = s-1; dt = 0
                    T0.remove([0,-e])
                    T0.remove([0,-e])
                    DT = tuple(sorted(T0))
                else:
                    dk = t+1; ds = 0; dt = 0
                    T0.remove([0,-e])
                    T0.remove([1,1])
                    for i in range(T0.count([0,-e])):
                        T0.remove([0,-e])
                        T0.append([0,1])
                    DT = tuple(sorted(T0))
        if T.count([1,-1]) > 0:
            if (T.count([0,e])-s) % 2 == 1: 
                if t % 2 == 0:
                    dk = t; ds = s; dt = 0 
                    DT = tuple(sorted(T0))
                else: 
                    dk = t; ds = s+1; dt = 0 
                    T0.append([0,e])
                    T0.remove([1,-1])
                    DT = tuple(sorted(T0))
            else:
                if t % 2 == 1:
                    dk = t+1; ds = s; dt = 0
                    T0.remove([0,e])
                    T0.remove([1,-1])
                    DT = tuple(sorted(T0)) 
                elif s > 0:
                    dk = t+1; ds = s-1; dt = 0
                    T0.remove([0,e])
                    T0.remove([0,e])
                    DT = tuple(sorted(T0))
                else:
                    dk = t+1; ds = 0; dt = 0
                    T0.remove([0,e])
                    T0.remove([1,-1])
                    for i in range(T0.count([0,e])):
                        T0.remove([0,e])
                        T0.append([0,-1])
                    DT = tuple(sorted(T0))
        if T.count([1,1])+T.count([1,-1]) == 0:    
            if T.count([0,1]) > 0 and (T.count([0,1])-s) % 2 == 1:          
                if t == 0:                       
                    dk = 0; ds = s; dt = 0       
                    DT = tuple(sorted(T0))          
                elif t % 2 == 0:         
                    dk = t-1; ds = s+1; dt = 0     
                    T0.append([0,1])             
                    T0.append([0,1])             
                    DT = tuple(sorted(T0))       
                else:             
                    dk = t-1; ds = s; dt = 1    
                    DT = tuple(sorted(T0))      
            if T.count([0,-1]) > 0 and (T.count([0,-1])-s) % 2 == 1:        
                if t == 0:                      
                    dk = 0; ds = s; dt = 0      
                    DT = tuple(sorted(T0))      
                elif t % 2 == 0:        
                    dk = t-1; ds = s+1; dt = 0    
                    T0.append([0,-1])           
                    T0.append([0,-1])           
                    DT = tuple(sorted(T0))      
                else:        
                    dk = t-1; ds = s; dt = 1            
                    DT = tuple(sorted(T0))                
            if T.count([0,1]) > 0 and (T.count([0,1])-s) % 2 == 0:
                if t % 2 == 1 and T.count([0,1]) > s and s == 0: 
                    dk = t; ds = s; dt = 0
                    for i in range(T0.count([0,1])):
                        T0.remove([0,1])
                        T0.append([0,-1])
                    DT = tuple(sorted(T0))
                elif t % 2 == 1 and T.count([0,1]) > s and s > 0:
                    dk = t; ds = s-1; dt = 1
                    T0.remove([0,1])
                    T0.remove([0,1])
                    DT = tuple(sorted(T0))
                else: 
                    dk = t; ds = s; dt = 0
                    DT = tuple(sorted(T0))                    
            if T.count([0,-1]) > 0 and (T.count([0,-1])-s) % 2 == 0:
                if t % 2 == 1 and T.count([0,-1]) > s and s == 0: 
                    dk = t; ds = s; dt = 0
                    for i in range(T0.count([0,-1])):
                        T0.remove([0,-1])
                        T0.append([0,1])
                    DT = tuple(sorted(T0))
                elif t % 2 == 1 and T.count([0,-1]) > s and s > 0:
                    dk = t; ds = s-1; dt = 1
                    T0.remove([0,-1])
                    T0.remove([0,-1])
                    DT = tuple(sorted(T0))
                else: 
                    dk = t; ds = s; dt = 0
                    DT = tuple(sorted(T0)) 
            if T.count([0,1])+T.count([0,-1]) == 0:
                dk = t; ds = 0; dt = 0             
                DT = tuple(sorted(T0))
        return dk, ds, dt, DT   


###########################################################
# The [0,1]-socles for G: A special case.
###########################################################

def S01A(k,s,t,T): 
    """
    [0,1]-socle: A special case. 
    Input: (k,s,t,T)
        k (non-negative integer)
        s (non-negative integer)
        t (non-negative integer)
        T (tempered L-parameter)
    Output (k, ds, dt, dT)
        ds (non-negative integer)
        dt (non-negative integer)                
        dT (tempered L-parameter)
            soc( Z[0,1]**{k} \rtimes L((||**1)**{s} \times [0,-1]**{t} \times T) )
            = L((||**1)**{ds} \times [0,-1]**{dt} \times dT)
    """
    n = []
    for j in range(s):                                        
        n.append([-1,-1])                                                        
    for j in range(t):                                        
        n.append([0,-1])
    m = tuple(n)                                      
    if D(1,m,T)[0] > 0:                                       
        print("L([-1,-1]**s,[0,-1]**t; T) has a nontrivial 1-derivative.")
    elif k == 0:
        return s,t,T
    else:
        T0 = list(T)    
        if k % 2 == 0:
            if t == 1:
                ds = s; dt = k+1; DT = tuple(sorted(T0))
            elif (T.count([0,1])+T.count([0,-1])-s) % 2 == 0:
                ds = s; dt = k; DT = tuple(sorted(T0))
            elif T.count([1,1])+T.count([1,-1]) > 0:
                ds = s; dt = k; DT = tuple(sorted(T0))
            else: 
                ds = s; dt = k-1
                if T.count([0,1]) > 0:
                    T0.append([0,1])
                    T0.append([1,(-1)**k])
                    DT = tuple(sorted(T0))
                if T.count([0,-1]) > 0:
                    T0.append([0,-1])
                    T0.append([1,-(-1)**k])
                    DT = tuple(sorted(T0))
        if k % 2 == 1:
            if t == 1:
                ds = s+1; dt = k
                if T.count([0,1]) > 0:
                    T0.append([0,1])
                    T0.append([0,1])
                    DT = tuple(sorted(T0))
                if T.count([0,-1]) > 0:
                    T0.append([0,-1])
                    T0.append([0,-1])
                    DT = tuple(sorted(T0))
            elif T.count([0,1])+T.count([0,-1]) == s:
                ds = s; dt = k; DT = T
            elif (T.count([0,1])+T.count([0,-1]) - s) % 2 == 0: 
                if s == 0:
                    ds = 0; dt = k; 
                    for i in range(T.count([0,1])):
                        T0.remove([0,1])
                        T0.append([0,-1])
                    for i in range(T.count([0,-1])):
                        T0.remove([0,-1])
                        T0.append([0,1])
                    DT = tuple(sorted(T0))
                else:
                    ds = s-1; dt = k+1
                    if T.count([0,1]) > 0:
                        T0.remove([0,1])
                        T0.remove([0,1])
                        DT = tuple(sorted(T0))
                    if T.count([0,-1]) > 0:
                        T0.remove([0,-1])
                        T0.remove([0,-1])
                        DT = tuple(sorted(T0))
            elif s == 0:
                if T.count([0,1]) > 0:
                    if T.count([1,1])+T.count([1,-1]):
                        ds = 1; dt = k-1
                        T0.append([0,1])
                        T0.append([0,1])
                    else:
                        ds = 0; dt = k-1
                        T0.append([0,1])
                        for i in range(T0.count([0,1])):
                            T0.remove([0,1])
                            T0.append([0,-1])
                        T0.append([1,-(-1)**k])
                    DT = tuple(sorted(T0))
                if T.count([0,-1]) > 0:
                    if T.count([1,1])+T.count([1,-1]):
                        ds = 1; dt = k-1
                        T0.append([0,-1])
                        T0.append([0,-1])
                    else:                    
                        ds = 0; dt = k-1
                        T0.append([0,-1])
                        for i in range(T0.count([0,-1])):
                            T0.remove([0,-1])
                            T0.append([0,1])
                        T0.append([1,(-1)**k])
                    DT = tuple(sorted(T0))                
            elif T.count([1,1])+T.count([1,-1]) > 0:
                ds = s+1; dt = k-1
                if T.count([0,1]) > 0:
                    T0.append([0,1])
                    T0.append([0,1])
                    DT = tuple(sorted(T0))
                if T.count([0,-1]) > 0:
                    T0.append([0,-1])
                    T0.append([0,-1])
                    DT = tuple(sorted(T0))
            else:                
                ds = s-1; dt = k
                if T.count([0,1]) > 0:
                    T0.remove([0,1])
                    T0.append([1,-(-1)**k])
                    DT = tuple(sorted(T0))
                if T.count([0,-1]) > 0:
                    T0.remove([0,-1])
                    T0.append([1,(-1)**k])
                    DT = tuple(sorted(T0))                
        return ds,dt,DT


###########################################################
# The [0,1]-derivatives for G: Reduction.
###########################################################

def D01(m,T):
    """
    Highest [0,1]-derivative. 
    Input: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (k,dm,dT)
        k (non-negative integer)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            D_{[0,1]}**{max}(m,T) = D_{x}**{(k)}(m,T) = (dm,dT)
    """
    if D(1,m,T)[0] > 0:
        print("L(m; T) has a nontrivial 1-derivative.")       
    else: 
        t = m.count([0,-1])
        s = m.count([-1,-1])
        n = list(filter(lambda z: z[1] != -1, m))
        n1 = []
        for i in range(s):
            n1.append([-1,-1])
        for i in range(t):
            n1.append([0,-1])
        l, m0, T0 = D(1,n1,T)
        kA, sA, tA, TA = D01A(m0.count([-1,-1]), m0.count([0,-1]),T0)
        for i in range(kA):
            n.append([0,0])
        for i in range(kA+l):
            n.append([1,1])
        k0,n = LD(0,n)
        k1,n = LD(1,n)
        for i in range(k0-k1):
            n = LS(0,n)
        n = list(n)
        k2 = n.count([0,0])
        k3 = n.count([1,1])
        l2 = k3-k2
        n = list(filter(lambda z: z != [0,0] and z != [1,1], n))
        s2,t2,T2 = S01A(k2,sA,tA,TA)
        n2 = []
        for i in range(s2):
            n2.append([-1,-1])
        for i in range(t2):
            n2.append([0,-1])
        for i in range(l2):
            n2,T2 = S(1,n2,T2)
        for i in range(n2.count([-1,-1])):
            n.append([-1,-1])
        for i in range(n2.count([0,-1])):
            n.append([0,-1])        
        dm = tuple(sorted(n))
        return k1, dm, T2


###########################################################
# The [0,1]-socles for G: Reduction.
###########################################################

def S01(k,m,T): 
    """
    [0,1]-socle.
    Input: (k,m,T)
        k (nonzero real number)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (dm,dT)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            soc( Z[0,1]**{k} \rtimes (m,T) ) = (dm,dT)
    """
    if D(1,m,T)[0] > 0:
        print("L(m; T) has a nontrivial 1-derivative.")       
    else:
        t = m.count([0,-1])
        s = m.count([-1,-1]) 
        n = list(filter(lambda z: z[1] != -1, m))
        n1 = []
        for i in range(s):
            n1.append([-1,-1])
        for i in range(t):
            n1.append([0,-1])
        l, m0, T0 = D(1,n1,T)
        kA, sA, tA, TA = D01A(m0.count([-1,-1]), m0.count([0,-1]),T0)
        for i in range(kA):
            n.append([0,0])
        for i in range(kA+l):
            n.append([1,1])
        k0,n = LD(0,n)
        for i in range(k):
            n = LS(1,n)
        for i in range(k+k0):
            n = LS(0,n)   
        n = list(n)
        k2 = n.count([0,0])
        k3 = n.count([1,1])
        l2 = k3-k2
        n = list(filter(lambda z: z != [0,0] and z != [1,1], n))
        s2,t2,T2 = S01A(k2,sA,tA,TA)
        n2 = []
        for i in range(s2):
            n2.append([-1,-1])
        for i in range(t2):
            n2.append([0,-1])
        for i in range(l2):
            n2,T2 = S(1,n2,T2)
        for i in range(n2.count([-1,-1])):
            n.append([-1,-1])
        for i in range(n2.count([0,-1])):
            n.append([0,-1])        
        dm = tuple(sorted(n))
        return dm, T2

###########################################################
# The Aubert duals for G.
###########################################################

def AD(m,T): 
    """
    Aubert dual.
    Input: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (dm,dT)
        dm (multi-segment, increasing negative exponents)
        dT (tempered L-parameter)
            Aubert dual of (m,T) = (dm,dT)
    """       
    if len(m) > 0 and len(T) > 0: 
        x0 = max([T[j][0] for j in range(len(T))])                                
        x1 = min([m[i][0] for i in range(len(m))])  
        x2 = max([-m[i][1] for i in range(len(m))]) 
        x = max(x0,x1,x2)                           
    if len(m) > 0 and len(T) == 0: 
        x1 = min([m[i][0] for i in range(len(m))])  
        x2 = max([-m[i][1] for i in range(len(m))]) 
        x = max(x1,x2)                                       
    if len(m) == 0 and len(T) > 0: 
        x0 = max([T[j][0] for j in range(len(T))])                                
        x = x0
    if len(m) == 0 and len(T) == 0: 
        return m, T                           
    while x > 0:                              
        if D(x,m,T)[0] > 0:                   
            break                             
        else:                                 
            x = x-1                           
    if x > 0:                                 
        k = D(x,m,T)[0]                       
        m1, T1 = AD(D(x,m,T)[1],D(x,m,T)[2])  
        for i in range(k):                    
            m1, T1 = S(-x,m1,T1)              
        return m1, T1                         
    elif len(m) > 0 and x1 < 0:               
        k = D(x1,m,T)[0]                      
        m1, T1 = AD(D(x1,m,T)[1],D(x1,m,T)[2])
        for i in range(k):           
            m1, T1 = S(-x1,m1,T1)               
        return m1, T1                         
    elif len(m) > 0:                          
        k, n, T = D00(m,T)
        m1, T1 = AD(n,T)                 
        m1, T1 = S01(k,m1,T1)                 
        return m1, T1                
    else:                 
        if (T.count([0,1])+T.count([0,-1])) % 2 == 1 or (2*T[0][0]) %2 == 1:
            return (), T
        else: 
            x0 = max([T[j][0] for j in range(len(T))])
            T0 = list(sorted(T)); m1 = ()
            if x0 > 0: 
                m1 = ([0,-x0],)
                T0.remove(T0[0])
                T0.remove(T0[len(T0)-1])
            for i in range(len(T0)):
                T0[i] = [T0[i][0],-T0[i][1]]
            T1 = tuple(T0)
            return m1, T1


###########################################################
# The multi-segment representations for G.
###########################################################

def rep(E):
    """
    Multi-segment representation.
    Input: E
        E (extended multi-segment)
    Output: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
            pi(E) = (m,T)
    """
    E0 = list(E)
    if len(E) == 0:
        return (),()
    else:
        b0 = min(E[i][0][1] for i in range(len(E)))
        if b0 < 0: 
            E0 = [ ([s[0][0]-int(b0-1),s[0][1]-int(b0-1)],s[1],s[2]) for s in E0]
            m, T = rep(E0)
            m0 = [ [z[0]+int(b0-1),z[1]-int(b0-1)] for z in m]
            T0 = [ [z[0]+int(b0-1),z[1]] for z in T]
        else:
            t = []
            for i in range(len(E0)):
                t.append(0)
            for i in range(len(E0)):
                if i > 0 and E0[i][0][1] < E0[i-1][0][0]+1:
                    t[i] = E0[i-1][0][0]-E0[i][0][1]+1
                    E0[i] = [[E0[i][0][0]+t[i],E0[i][0][1]+t[i]],E0[i][1],E0[i][2]]
            m0 = []; T0 = []
            for i in range(len(E0)):
                A = E0[i][0][0]; B = E0[i][0][1]; l = 0
                while l < E0[i][1]:
                    m0.append([B+l,-A+l])
                    l = l+1
            for i in range(len(E0)):
                A = E0[i][0][0]-E0[i][1]; B = E0[i][0][1]+E0[i][1]
                l = 0; e = E0[i][2]
                while l <= A-B:
                    T0.append([B+l,e])
                    l = l+1
                    e = -e
            for i in range(len(E0)):
                A = E0[i][0][0]; B = E0[i][0][1]; s = 0
                while s < t[i]:
                    l = 0
                    while l <= A-B:
                        m0,T0 = D1(B+l-s,m0,T0)
                        l = l+1
                    s = s+1
        for j in range(len(m0)):
            if m0[j][0] < m0[j][1]-1:
                m0 = (); T0 = ()
                break
        for j in range(len(T0)):
            if T0[j][0] < -1/2 or T0[j] == [-1/2,-1]:
                m0 = (); T0 = ()
                break
        dm = tuple(sorted(filter(lambda z: z[0] >= z[1], m0)))
        dT = tuple(sorted(filter(lambda z: z != [-1/2,+1], T0)))
        return dm,dT


###########################################################
# Xu’s algorithm for rep(E) != 0.
###########################################################

def nec(E):
    r"""
    Necessary condition for \pi(E) != 0.
    Input: E
        E (extended multi-segment)
    Output: True or False 
    """
    e = 1
    for i in range(len(E)-1):
        A1 = E[i][0][0]; B1 = E[i][0][1]; l1 = E[i][1]; e1 = E[i][2]
        A2 = E[i+1][0][0]; B2 = E[i+1][0][1]; l2 = E[i+1][1]; e2 = E[i+1][2]        
        if A2 >= A1 and B2 >= B1:
            if e2*e1 == (-1)**(A1-B1):
                if A2-l2 < A1-l1 or B2+l2 < B1+l1:
                    e = 0
                    break
            else:
                if B2+l2 <= A1-l1:
                    e = 0
                    break
        if A2 >= A1 and B2 <= B1:
            if e2*e1 == (-1)**(A1-B1):
                if l2-l1 < 0 or l2-l1 > (A2-B2)-(A1-B1):
                    e = 0
                    break
            else:
                if l2+l1 < A1-B1+1:
                    e = 0
                    break
        if A2 <= A1 and B2 >= B1:
            if e2*e1 == (-1)**(A1-B1):
                if l1-l2 < 0 or l1-l2 > (A1-B1)-(A2-B2):
                    e = 0
                    break
            else:
                if l1+l2 < A2-B2+1:
                    e = 0
                    break
    return e == 1


def change(E,i):
    """
    Changing admissible orders.
    Input: (E,i)
        E (extended multi-segment)
        i (integer with 0 <= i <= len(E)-2)
    Output: E1
        E1 (extended multi-segment)
            change order so that i > i+1 if the resulting order is admissible.
    """
    E0 = list(E)
    A1 = E[i][0][0]; B1 = E[i][0][1]; l1 = E[i][1]; e1 = E[i][2]
    A2 = E[i+1][0][0]; B2 = E[i+1][0][1]; l2 = E[i+1][1]; e2 = E[i+1][2]   
    if A2 >= A1 and B2 <= B1:
        if e2*e1 == (-1)**(A1-B1+1):
            E0[i] = ([A2,B2],l2+2*l1-(A1-B1+1),e1)
            E0[i+1] = ([A1,B1],l1,(-1)**(A2-B2)*e1)
        else: 
            if l2-2*l1 < (A2-B2)/2-(A1-B1):
                E0[i] = ([A2,B2],l2-2*l1+(A1-B1+1),-e1)
                E0[i+1] = ([A1,B1],l1,(-1)**(A2-B2)*e1)
            else: 
                E0[i] = ([A2,B2],-l2+2*l1+(A2-B2)-(A1-B1),e1)
                E0[i+1] = ([A1,B1],l1,(-1)**(A2-B2)*e1)
    if A2 <= A1 and B2 >= B1:
        if e2*e1 == (-1)**(A1-B1+1):
            E0[i] = ([A2,B2],l2,-e1)
            E0[i+1] = ([A1,B1],l1+2*l2-(A2-B2+1),(-1)**(A2-B2+1)*e1)
        else: 
            if l1-2*l2 < (A1-B1)/2-(A2-B2):
                E0[i] = ([A2,B2],l2,e1)
                E0[i+1] = ([A1,B1],l1-2*l2+(A2-B2+1),(-1)**(A2-B2+1)*e1)
            else: 
                E0[i] = ([A2,B2],l2,e1)
                E0[i+1] = ([A1,B1],-l1+2*l2+(A1-B1)-(A2-B2),(-1)**(A2-B2)*e1)
    return tuple(E0)

def orders(E): 
    """
    The set of all extended multi-segments 
    which are given from $E$ by changing admissible orders.
    Input: E
        E (extended multi-segment)
    Output: C_E 
        C_E (a set of extended multi-segments)
            C_E = {E' | E' is given from E by changing admissible orders}.
    """
    if len(E) <= 1:
        return([E])
    else: 
        result = []; L = []; top = []
        for i in range(len(E)):
            l = 0
            for j in range(len(E)):
                if E[j][0][0] < E[i][0][0] and E[j][0][1] < E[i][0][1]:
                    l = 1
                    break
            if l==0:
                k = 0
                for j in top: 
                    if E[i][0] == E[j][0]:
                        k=1
                        break
                if k == 0:
                    top.append(i)
        for i in top:
            E0 = list(E)
            if i == 0:
                L.append(E0)
            else:
                for j in range(i):
                    E0 = change(E0, i-j-1)
                L.append(E0)
        for E0 in L:
            E1 = list(E0)
            del E1[0]
            E1 = tuple(E1)
            for E2 in orders(E1):
                result.append(tuple([E0[0]] + list(E2)))
        return(result)

def nonzero(E,e):         
    s = 1
    for i in range(len(E)): 
        b = E[i][0][0]-E[i][0][1]+1
        s = (-1)**(int(b/2) + E[i][1]) * E[i][2]**b * s
    if s != e: 
            return False
    for i in range(len(E)): 
        if E[i][0][1]+E[i][1] <= -1: 
            return False
    for i in range(len(E)):
        d = 0
        for k in range(i): 
            d = d + (E[k][0][0]+E[k][0][1]+1)
        if (E[i][0][1] + E[i][1]) == -1/2 and E[i][2] != (-1)**d: 
            return False
    for E0 in orders(E):
        if not nec(E0):
            return False
    return True


###########################################################
# Deformation.
###########################################################

def deform(E,i): 
    """
    Deformation.
    Input: (E,i)
        E (extended multi-segment)
        i (integer with 0 <= i <= len(E))
    Output: E0
        E0 (extended multi-segment)
    """
    E0 = list(E)
    A1 = E[i][0][0]; B1 = E[i][0][1]; l1 = E[i][1]; e1 = E[i][2]
    if l1 == 0:
        if i < len(E)-1 and E[i+1][0][1] == A1+1 and E[i+1][1] == 0 and e1*E[i+1][2] == -(-1)**(A1-B1):
            del E0[i]
            del E0[i]
            E0.insert(i, ([E[i+1][0][0],B1],0,e1))
        else: 
            J1 = list(filter(lambda j: E[j][0][1] > B1, range(i)))
            if J1 == []:
                A0 = B1
            else: 
                A0 = max([E[j][0][0] for j in J1])
            J2 = list(filter(lambda j: E[j][0][0] < A1, range(i+1,len(E))))
            if J2 == []:
                B0 = A1
            else: 
                B0 = min([E[j][0][1] for j in J2])
            if A0 < B0:
                del E0[i]
                E0.insert(i, ([A0,B1],0,e1))
                E0.insert(i+1, ([A1,B0],0,(-1)**(B0-B1)*e1))
                for k in range(int(B0-A0-1)):
                    E0.insert(i+1,([B0-1-k,B0-1-k],0,(-1)**(B0-B1+1+k)*e1))
                if E0[i][0] == [-1/2,-1/2]:
                    for k in range(i):
                        E0 = list(change(E0,i-1-k))
                    del E0[0]
    elif i < len(E)-1:
        A1 = E[i][0][0]; B1 = E[i][0][1]; l1 = E[i][1]; e1 = E[i][2]
        A2 = E[i+1][0][0]; B2 = E[i+1][0][1]; l2 = E[i+1][1]; e2 = E[i+1][2]  
        if B1 < B2 and A2 < A1:
            J1 = list(filter(lambda j: B1 < E[j][0][1] and A2 < E[j][0][0], range(i)))
            J2 = list(filter(lambda j: E[j][0][1] < B2 and E[j][0][0] < A1, range(i+2,len(E))))
            if J1 == [] and J2 == []:
                if 2*l2 == A2-B2+1:
                    if l1 == l2:
                        E0[i] = ([A2,B1],l1,e1)
                        E0[i+1] = ([A1,B2],l1,-(-1)**(A2-B1)*e1)
                        if A2+B1 < 0:
                            for k in range(i):
                                E0 = list(change(E0,i-1-k))
                            del E0[0]
                elif 2*l1 == A1-B1+1:
                    if l1 == (A2-B2+1)-l2:
                        E0[i] = ([A2,B1],l2+(B2-B1),(-1)**(A1-B1)*e2)
                        E0[i+1] = ([A1,B2],l2,(-1)**(A2-A1)*e2)
                        if A2+B1 < 0:
                            for k in range(i):
                                E0 = list(change(E0,i-1-k))
                            del E0[0]
                elif e1*e2 == (-1)**(A1-B1) and l1 == l2:
                    if 2*l1+A1-A2 <= A2-B2+1:
                        E0[i] = ([A2,B1],l1,e1)
                        E0[i+1] = ([A1,B2],l1+(A1-A2),(-1)**(A1-A2)*e2)
                    else: 
                        E0[i] = ([A2,B1],l1,e1)
                        E0[i+1] = ([A1,B2],-l1+(A2-B2+1),-(-1)**(A1-A2)*e2)
                    if A2+B1 < 0:
                        for k in range(i):
                            E0 = list(change(E0,i-1-k))
                        del E0[0]
                elif e1*e2 == -(-1)**(A1-B1) and l1 == (A2-B2+1)-l2:
                    if l1+B1 <= l2+B2:
                        E0[i] = ([A2,B1],l1,e1)
                        E0[i+1] = ([A1,B2],l2,-(-1)**(A2-B1)*e1)
                    else:
                        E0[i] = ([A2,B1],l2+(B2-B1),-e1)
                        E0[i+1] = ([A1,B2],l2,-(-1)**(A2-B1)*e1)
                    if A2+B1 < 0:
                        for k in range(i):
                            E0 = list(change(E0,i-1-k))
                        del E0[0]
        if B1 < B2 and A1 < A2: 
            if A2-l2 == A1-l1 and e2*e1 == (-1)**(A1-B1):
                E0[i] = ([A2,B1],l1,e1)
                E0[i+1] = ([A1,B2],l2-(A2-A1),(-1)**(A2-A1)*e2)
                if B2 == A1+1:
                   del E0[i+1]
            if B2+l2 == B1+l1 and e2*e1 == (-1)**(A1-B1):
                if (A1-B1+1)-2*l1 >= A2-A1:
                    E0[i] = ([A2,B1],l1+(A2-A1),e1)
                    E0[i+1] = ([A1,B2],l2,(-1)**(A2-A1)*e2)
                else: 
                    E0[i] = ([A2,B1],(A1-B1+1)-l1,-e1)
                    E0[i+1] = ([A1,B2],l2,(-1)**(A2-A1)*e2)
                if B2 == A1+1:
                   del E0[i+1]
            if B2+l2 == A1-l1+1 and e1*e2 == -(-1)**(A1-B1):
                if l2 <= l1:
                    E0[i] = ([A2,B1],l1,e1)
                    E0[i+1] = ([A1,B2],l2,(-1)**(A2-A1)*e2)
                else: 
                    E0[i] = ([A2,B1],l1,e1)
                    E0[i+1] = ([A1,B2],l1,-(-1)**(A2-A1)*e2)
                if B2 == A1+1:
                   del E0[i+1]
    return tuple(E0)

###########################################################
# A-packets for G.
###########################################################

def A_packet(P,e):
    """
    A-packet associated to P.
    Input: (P,e)
        P (A-parameter)
        e (sign +1 or -1)
    Output: Pi
        Pi (A-packet)
    """
    Pi = []; S = []; T = []; E = []; l = []
    for i in range(len(P)):
        s = (P[i][0]-P[i][1], [(P[i][0]+P[i][1])/2-1, (P[i][0]-P[i][1])/2])
        S.append(s)
    S = sorted(S)
    t = 0
    for i in range(len(S)):
        if i == 0:
            T.append([S[i][1],1])
            l.append(1)
        elif S[i] != S[i-1]:
            T.append([S[i][1],1])
            t = t+1
            l.append(1)
        elif S[i] == S[i-1]:
            T[t] = [T[t][0], T[t][1]+1]
    for i in range(len(T)):
        for j in range(T[i][1]):
            E.append( (T[i][0], 0, (-1)**(j*(T[i][0][0]-T[i][0][1]))) )
    E = tuple(E)
    Pi.append(E)
    t = 0
    for i in range(len(T)): 
        for j in range(len(Pi)): 
            F = list(Pi[j])
            while l[i] <= (T[i][0][0]-T[i][0][1]+1)/2: 
                if 2*l[i] < T[i][0][0]-T[i][0][1]+1:
                    for k in range(T[i][1]): 
                        F[t+k] = (F[t+k][0], l[i], F[t+k][2])
                if 2*l[i] == T[i][0][0]-T[i][0][1]+1:
                    for k in range(T[i][1]): 
                        F[t+k] = (F[t+k][0], l[i], 1)
                Pi.append(tuple(F))
                l[i] = l[i]+1
            l[i] = 1
        t = t+T[i][1]
    t = 0
    for i in range(len(T)): 
        for j in range(len(Pi)): 
            F = list(Pi[j])
            if F[t][2] == +1 and F[t][1] < (F[t][0][0]-F[t][0][1]+1)/2:
                for k in range(T[i][1]):
                    F[t+k] = (F[t+k][0], F[t+k][1], -F[t+k][2])
                Pi.append(tuple(F))
        t = t+T[i][1]
    Pi = list(filter(lambda F: nonzero(F,e), Pi))
    return Pi


def L_packet(phi,e):
    """
    L-packet associated to phi
    Input: (phi,e)
        phi (tempered L-parameter)
        e (sign +1 or -1)
    Output: Pi
        Pi (tempered L-packet)
    """
    return [rep(E)[1] for E in A_packet(phi,e)] 


###########################################################
# The symbols of extended multi-segments.
###########################################################

def symbol(E):
    """
    Symbol.
    Input: E
        E (extended multi-segment)
    Output: printing symbol(E)
    """
    n = len(E)
    A0 = max([ E[i][0][0] for i in range(n) ])
    B0 = min([ E[i][0][1] for i in range(n) ])
    X = np.array([B0+i for i in range(A0-B0+1)])  
    for i in range(n):
        v = []
        A = E[i][0][0]
        B = E[i][0][1]
        l = E[i][1]
        e = E[i][2]
        for j in range(A0-B0+1):
            if B0+j < B or B0+j > A:
                if B0+j < 0 and (2*B0) % 2 == 0:
                    v.append("  ")
                if B0+j >= 0 and (2*B0) % 2 == 0:
                    v.append(" ")
                if B0+j < 0 and (2*B0) % 2 == 1:
                    v.append("    ")
                if B0+j >= 0 and (2*B0) % 2 == 1:
                    v.append("   ")
            elif B0+j-B < l:
                if B0+j < 0 and (2*B0) % 2 == 0:
                    v.append(" <")
                if B0+j >= 0 and (2*B0) % 2 == 0: 
                    v.append("<")
                if B0+j < 0 and (2*B0) % 2 == 1:
                    v.append("  < ")
                if B0+j >= 0 and (2*B0) % 2 == 1: 
                    v.append(" < ")
            elif A-(B0+j) < l:
                if B0+j < 0 and (2*B0) % 2 == 0:
                    v.append(" >")
                if B0+j >= 0 and (2*B0) % 2 == 0:
                    v.append(">")
                if B0+j < 0 and (2*B0) % 2 == 1:
                    v.append("  > ")
                if B0+j >= 0 and (2*B0) % 2 == 1:
                    v.append(" > ")
            elif (-1)**(B0+j-l-B) == e:
                if B0+j < 0 and (2*B0) % 2 == 0:
                    v.append(" +")
                if B0+j >= 0 and (2*B0) % 2 == 0:
                    v.append("+")
                if B0+j < 0 and (2*B0) % 2 == 1:
                    v.append("  + ")
                if B0+j >= 0 and (2*B0) % 2 == 1:
                    v.append(" + ")
            else:
                if B0+j < 0 and (2*B0) % 2 == 0:
                    v.append(" -")
                if B0+j >= 0 and (2*B0) % 2 == 0:
                    v.append("-")
                if B0+j < 0 and (2*B0) % 2 == 1:
                    v.append("  - ")
                if B0+j >= 0 and (2*B0) % 2 == 1:
                    v.append(" - ")
        X = np.append(X,v)
    X = X.reshape(n+1,A0-B0+1)
    print(X)
#    return X

###########################################################
# The associated A-parameter for E.
###########################################################

def dim(psi):
    """
    Dimension of an A-parameter
    Input: psi
        psi (A-parameter)
    Output: d
        d (non-negative integer)
            d = dim(psi)
    """
    d = 0
    for a in psi:
        d = d + a[0]*a[1]
    return d

def par(E):
    """
    Parameter of extended multi-segments.
    Input: E
        E (extended multi-segment)
    Output: psi
        psi (A-parameter)
            psi = psi_E
    """
    psi = []
    for S in E: 
        a = S[0][0]+S[0][1]+1
        b = S[0][0]-S[0][1]+1
        psi.append((a,b))
    return tuple(psi)

###########################################################
# Character for rep(E).
###########################################################

def char(E): 
    """
    character associated to E
    Input: E
        E (extended multi-segment)
    Output: (C, psi)
        C (list of +1 or -1)
        psi (A-parameter associated to E)
    """
    C = []; z = []
    for i in range(len(E)):
        C.append(1); z.append(0) 
    for i in range(len(E)):
        a1 = E[i][0][0]+E[i][0][1]+1; b1 = E[i][0][0]-E[i][0][1]+1
        for j in range(len(E)):
            a2 = E[j][0][0]+E[j][0][1]+1; b2 = E[j][0][0]-E[j][0][1]+1
            if b1%2 == 0 and b2%2 == 1:
                if j<i and a2 > a1 and b1 > b2:
                    z[i] = z[i]+1
                if j>i and a1 > a2 and b1 > b2:
                    z[i] = z[i]+1
            if b1%2 == 1 and b2%2 == 0:
                if j<i and a2 > a1 and b2 > b1:
                    z[i] = z[i]+1
                if j>i and a1 > a2 and b2 > b1:
                    z[i] = z[i]+1
        C[i] = (-1)**(z[i]+int(b1/2)+E[i][1])*E[i][2]**(b1)
    return C, par(E)


###########################################################
# Dual of E.
###########################################################

def hat(E):
    """
    Dual of an extended multi-segment. 
    Input: E
        E (extended multi-segment)
    Output: E0
        E0 (extended multi-segment)
    """
    E0 = list(E)
    F = []
    for i in range(len(E0)): 
        if (2*E0[i][0][1]) %2 == 0:
            b = 0
            for j in range(len(E0)):
                b = b + (E0[j][0][0]-E0[j][0][1]+1)
            F.insert(0,([E0[i][0][0],-E0[i][0][1]],E0[i][1]+E0[i][0][1],(-1)**(b-E0[i][0][0]+E0[i][0][1]-1)*E0[i][2])) 
        if (2*E0[i][0][1]) %2 == 1:
            a = 0; b = 0
            for j in range(len(E0)): 
                if j > i: 
                   a = a + (E0[j][0][0]+E0[j][0][1])
                if j < i:
                   b = b + (E0[j][0][0]-E0[j][0][1])
            if E[i][2] != (-1)**b or E0[i][1] == (E0[i][0][0]-E0[i][0][1]+1)/2: 
                F.insert(0,([E[i][0][0],-E[i][0][1]],E[i][1]+E[i][0][1]-1/2,(-1)**a))
            else:
                F.insert(0,([E[i][0][0],-E[i][0][1]],E[i][1]+E[i][0][1]+1/2,(-1)**(a+1)))
    return tuple(F)

def dual_rep(E):
    """
    Aubert dual of rep(E).
    Input: E
        E (extended multi-segment)
    Output: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    """
    return rep(hat(E))


###########################################################
# Strongly equivalence classes.
###########################################################

def reorder(E):
    """
    Reorder an extended multi-segment.
    Input: E
        E (extended multi-segment)
    Output: E0
        E0 (extended multi-segment)
    """
    E0 = list(E)
    a = 1
    while a == 1:
        a = 0
        for i in range(len(E0)-1):
            if E0[i][0][1] > E0[i+1][0][1] and E0[i][0][0] <= E0[i+1][0][0]:
                E0 = list(change(E0,i))
                a = 1
                break
    b = 1
    while b == 1:
        b = 0
        for i in range(len(E0)-1):
            if E0[i][0][1] == E0[i+1][0][1] and E0[i][0][0] > E0[i+1][0][0]:
                E0 = list(change(E0,i))
                b = 1
                break
    t=0
    for i in range(len(E0)):
        if E0[i-t][0][0]+E0[i-t][0][1] < 0:
            for k in range(i-t):
                E0 = list(change(E0,i-t-1-k))
            del E0[0]
            t += 1
    for i in range(len(E0)-1):
        if 2*E0[i][1] == E0[i][0][0]-E0[i][0][1]+1:
            E0[i] = (E0[i][0],E0[i][1],+1)
    return tuple(E0)


def CUIP(E,C):
    """
    Add extended multi-segments given from E by three operators to C.
    Input: (E,C)
        E (extended multi-segment)
        C (list)
    Output: C1
        C1 (list)
    """
    if E == (): 
        return C
    else: 
        E0 = tuple(E)
        e = (2*E[0][0][0] % 2)/2
        A = max([z[0][0] for z in E])
        m = len(E)
        for k in range(m):
            E1 = list(E0)
            E1 = deform(E1,k)
            E1 = reorder(E1)
            if E1 not in C:
                C.append(E1) 
        for l in range(int(-E0[0][0][1]+1-e),int(E0[0][0][0]+1-e)):
            E1 = list(E0)
            E1.insert(0,([l-1+e,-l-e],l,+1))
            E1 = deform(E1,0)
            E1 = reorder(E1)
            if E1 not in C:
                C.append(E1)
        for j in range(1,m):
            E1 = list(E0)
            for j1 in range(j):
                if E1 == change(E1,j-j1-1):
                    break
                else:
                    E1 = change(E1,j-j1-1)
                    E1 = tuple(E1)
                    if j1 < j-1:
                        E2 = list(E1)
                        E2 = deform(E2,j-j1-2)
                        E2 = reorder(E2)
                        if E2 not in C:
                            C.append(E2)
                    else:
                        for l in range(int(-E1[0][0][1]+1-e),int(E1[0][0][0]+1-e)):
                            E2 = list(E1)
                            E2.insert(0,([l-1+e,-l-e],l,+1))
                            E2 = deform(E2,0)
                            E2 = reorder(E2)
                            if E2 not in C:
                                C.append(E2)
        for i in range(m-1):
            E1 = list(E0)
            for i1 in range(m-i-2):
                if E1 == change(E1,i+i1):
                    break
                else:
                    E1 = change(E1,i+i1)
                    E1 = tuple(E1)
                    for j in range(i+i1+2,m):
                        E2 = list(E1)
                        for j1 in range(j-(i+i1+2)):
                            if E2 == change(E2,j-j1-1):
                                break
                            else:
                                E2 = change(E2,j-j1-1)
                        if j == i+i1+2 or j1 == j-(i+i1+2)-1:
                            E2 = deform(E2,i+i1+1)
                            E2 = reorder(E2)
                            if E2 not in C:
                                C.append(E2)
        for l in range(-1,int(A+1-e)):
            if l == -1 or (l == 0 and e == 0):
                pass
            else: 
                E1 = list(E0)
                E1.insert(0,([l-1+e,-l-e],l,+1))
                E1 = tuple(E1)
                for i1 in range(m-1):
                    if E1 == change(E1,i1):
                        break
                    else:
                        E1 = change(E1,i1)
                        E1 = tuple(E1)
                        for j in range(i1+2,m):
                            E2 = list(E1)
                            for j1 in range(j-(i1+2)):
                                if E2 == change(E2,j-j1-1):
                                    break
                                else:
                                    E2 = change(E2,j-j1-1)
                            if j == i1+2 or j1 == j-(i1+2)-1:
                                E2 = deform(E2,i1+1)
                                E2 = reorder(E2)
                                if E2 not in C:
                                    C.append(E2)
    return C

def eq_class(E):
    """
    Strongly equivalence class. 
    Input: E
        E (extended multi-segment)
    Output: C
        C (list)
            C = [ E1 : rep(E1) == rep(E) ]
    """
    C = [E]
    C1 = []
    while len(C)-len(C1) > 0:
        for F in C:
            if F not in C1:
                C = CUIP(F,C)
                C1.append(F)
                break
    return C


###########################################################
# Is_Arthur?
###########################################################

def Is_Arthur(m,T):
    """
    Determination of Arthur type representations.
    Input: (m,T)
        m (multi-segment, increasing negative exponents)
        T (tempered L-parameter)
    Output: (bool, n, E)
        bool (True or False)
        n (non-negative integer)
        E (extended multi-segment)
            bool = True <=> (m,T) is of Arthur type. 
            In this case, (m,T) = rep(E), n = len(eq_class(E)).
    """
    if m == ():
        E = tuple([ ([p[0],p[0]],0,p[1]) for p in T])
        return True, len(eq_class(E)), E
    else:
        if T == ():
            x1 = min([s[0] for s in m])  
            x2 = max([-s[1] for s in m]) 
            x = max(x1,x2)  
        else:
            x0 = max([p[0] for p in T])                                
            x1 = min([s[0] for s in m])  
            x2 = max([-s[1] for s in m]) 
            x = max(x0,x1,x2)   
        while x > 0:
            k1,m1,T1 = D(x,m,T)
            if k1 > 0:
               break
            else: 
                x = x-1
        if x >= 1:
            V = [(x,k1)]
            i = 1
            while k1 > 0:
                k1,m1,T1 = D(x+i,m1,T1)
                if k1 > 0:
                    V.append((x+i,k1))
                    i += 1
            r = len(V)
            V.append((x+r,0))
            frag, _, E0 = Is_Arthur(m1,T1)
            if frag:
                f = 0
                C = eq_class(E0)
                for F in C:
                    S = [s[0] for s in F]
                    if F != () and max([s[1] for s in S]) == x-1:
                        if f == 1:
                            break
                        else: 
                            for i in range(r):
                                if S.count([x+r-2-i, x-1]) < V[r-1-i][1]-V[r-i][1]: 
                                    f = 0
                                    break 
                                else: 
                                    E1 = list(F)
                                    f = 1
                if f == 0:
                    return False, 0, ()
                else:
                    for i in range(r):
                        for _ in range(V[r-1-i][1]-V[r-i][1]):
                            J = filter(lambda j: E1[j][0][1] == x-1 and E1[j][0][0] == x+r-2-i, range(len(E1)))
                            j0 = max(list(J))
                            E1[j0] = ([x+r-1-i,x],E1[j0][1],E1[j0][2])
                    E1 = reorder(E1)
                    E = tuple(E1)
                    return True, len(eq_class(E)), E
            else:
                return False, 0, ()
        elif x1 < 0:
            x = x1
            k1,m1,T1 = D(x,m,T)
            V = [(x,k1)]
            i = 1
            while k1 > 0:
                k1,m1,T1 = D(x-i,m1,T1)
                if k1 > 0:
                    V.append((x-i,k1))
                    i += 1
            r = len(V)
            V.append((x-r,0))
            frag, _, E0 = Is_Arthur(m1,T1)
            if frag:
                if x == -1/2 and r == 1:
                    f = 1
                    E1 = list(E0)
                else: 
                    f = 0
                    C = eq_class(E0)
                    for F in C:
                        S = [s[0] for s in F]
                        if F != () and min([s[1] for s in S]) == x+1:
                            if f == 1:
                                break
                            else: 
                                for i in range(r):
                                    if x == -1/2 and i == r-1:
                                        pass
                                    elif S.count([-x+r-2-i,x+1]) < V[r-1-i][1]-V[r-i][1]:
                                        f = 0
                                        break 
                                    else: 
                                        E1 = list(F)
                                        f = 1
                if f == 0:
                    return False, 0, ()
                else:
                    for i in range(r):
                        for _ in range(V[r-1-i][1]-V[r-i][1]):
                            if x == -1/2 and i == r-1: 
                                E1.insert(0,([1/2,-1/2],1,1))
                            else: 
                                J = filter(lambda j: E1[j][0][1] == x+1 and E1[j][0][0] == -x+r-2-i, range(len(E1)))
                                j0 = min(list(J))
                                E1[j0] = ([-x+r-1-i,x],E1[j0][1]+1,E1[j0][2])
                    E1 = reorder(E1)
                    E = tuple(E1)
                    return True, len(eq_class(E)), E
            else:
                return False, 0, ()
        else: 
            S = []
            for s in m:
                S.append(s[0])
                S.append(-s[1])
            for p in T:
                S.append(p[0])            
            A = max(S)-m[0][0]
            k = [S.count(i+m[0][0]) for i in range(int(A)+1)]
            k.append(0)
            f = 1
            for i in range(int(A)):
                if k[i] < k[i+1]:
                    f = 0
                    break
            if f == 1:
                if x == 0:
                    psi0 = []
                    for i in range(int(A)+1):
                        for _ in range(k[i]-k[i+1]):
                            psi0.append((i+1,i+1))
                    psi = tuple(psi0)
                    Pi = A_packet(psi,+1)
                    f1 = 0
                    for E in Pi:
                        if rep(E) == (m,T):
                            f1 = 1
                            break
                    if f1 == 1:
                        return True, len(eq_class(E)), E
                    else: 
                        return False, 0, ()
                else:
                    psi0 = []
                    for i in range(int(A)+1):
                        for _ in range(k[i]-k[i+1]):        
                            psi0.append((i+2,i+1))
                    psi = tuple(psi0)
                    Pi = A_packet(psi,+1)
                    f1 = 0
                    for E in Pi:
                        if rep(E) == (m,T):
                            f1 = 1
                            break
                    if f1 == 1:
                        return True, len(eq_class(E)), E
                    else: 
                        return False, 0, ()
            else: 
                return False, 0, ()




###########################################################
# decomposition of Speh times Arthur.
###########################################################

def decomp(ab,E):
    r"""
    Decomposition of Speh(a,b) rtimes rep(E).
    Input: (ab,E)
        ab (= (a,b), two positive integers)
        E (extended multi-segment)
    Output: Pi
        Pi (list of extended multi-segments)
            Speh(a,b) rtimes rep(E) = oplus_{E1 in Pi} rep(E1).
    """
    a = ab[0]; b = ab[1]
    A = (a+b)/2-1; B = (a-b)/2
    c = 0
    for c in range(len(E)):
        if E[c][0][1] <= B:
            c = c+1
        else:
            break
    e = 0
    for i in range(c):
        e = e + (E[i][0][0]+E[i][0][1]+1)
    Pi = []
    for l in range(b//2+1):
        if B+l >= 0 or (B+l == -1/2 and e % 2 == 0):
            E0 = list(E)
            E0.insert(c,([A,B],l,(-1)**(b-1)))
            E0.insert(c,([A,B],l,1))
            E0 = tuple(E0)
            Pi.append(E0)
    for l in range((b+1)//2):
        if B+l >= 0 or (B+l == -1/2 and e % 2 == 1):
            E0 = list(E)
            E0.insert(c,([A,B],l,(-1)**b))
            E0.insert(c,([A,B],l,-1))
            E0 = tuple(E0)
            Pi.append(E0) 
    s = 1
    for i in range(len(E)): 
        b = E[i][0][0]-E[i][0][1]+1
        s = (-1)**(int(b/2) + E[i][1]) * E[i][2]**b * s
    Pi = list(filter(lambda F: nonzero(F,s), Pi))
    return Pi


###########################################################
# socle of Speh times Arthur.
###########################################################

def soc(s,ab,E):
    r"""
    Socle of Speh(a,b)||**{s} rtimes rep(E).
    Input: (s,ab,E)
        s (real number)
        ab (= (a,b), two positive integers)
        E (extended multi-segment)
    Output: Pi
        Pi (list of (m,T))
            soc(Speh(a,b)||**{s} rtimes rep(E)) = oplus_{(m,T) in Pi} (m,T).
    """
    a = ab[0]; b = ab[1]
    Pi = []
    if s == 0:
        Pi1 = decomp((a,b),E)
        for E1 in Pi1:
            Pi.append(rep(E1))
        return Pi
    A = (a+b)/2-1; B = (a-b)/2
    if s > (a-1)/2:
        m,T = rep(E)
        C = []
        for i in range(a):
            for j in range(b):
                C.append(B+s-i+j)
        L = []
        i = 0
        while i < len(C):
            if C[i] == 0:
                k,m,T = D(1,m,T)
                L.append(k)
                k,m,T = D01(m,T)
                L.append(k)
                i = i+2
            else:
                k,m,T = D(C[i],m,T)
                L.append(k)
                i = i+1
        C = list(reversed(C)); L = list(reversed(L))
        i = 0
        while i < len(C):
            if i < len(C)-1 and C[i+1] == 0:
                m,T = S01(L[i]+1,m,T)
                for j in range(L[i+1]):
                    m,T = S(1,m,T)
                i = i+2
            else:
                for j in range(L[i]+1):
                    m,T = S(C[i],m,T)
                i = i+1
        Pi.append((m,T))
        return Pi
    if s < -(b-1)/2:
        m,T = rep(E)
        C = []
        for j in range(b):
            for i in range(a):
                C.append(B+s-i+j)
        L = []
        i = 0
        while i < len(C):
            if C[i] == 0:
                k,m,T = D(-1,m,T)
                L.append(k)
                k,m,T = D00(m,T)
                L.append(k)
                i = i+2
            else:
                k,m,T = D(C[i],m,T)
                L.append(k)
                i = i+1
        C = list(reversed(C)); L = list(reversed(L))
        i = 0
        while i < len(C):
            if i < len(C)-1 and C[i+1] == 0:
                for j in range(L[i]+1):
                    m,T = S00(m,T)
                for j in range(L[i+1]):
                    m,T = S(-1,m,T)
                i = i+2
            else:
                for j in range(L[i]+1):
                    m,T = S(C[i],m,T)
                i = i+1
        Pi.append((m,T))
        return Pi
    if s > 0 and s <= (a-1)/2:
        m,T = rep(E)
        C = []
        for i in range(int(2*s)):
            for j in range(b):
                C.append(B+s-i+j)
        L = []
        i = 0
        while i < len(C):
            if C[i] == 0:
                k,m,T = D(1,m,T)
                L.append(k)
                k,m,T = D01(m,T)
                L.append(k)
                i = i+2
            else:
                k,m,T = D(C[i],m,T)
                L.append(k)
                i = i+1
        Cr = list(reversed(C)); Lr = list(reversed(L))                
        Pi1 = decomp((int(a-2*s),b),E)
        for E1 in Pi1:
            m1,T1 = rep(E1)
            L1 = []
            i = 0
            while i < len(C):
                if C[i] == 0:
                    k1,m1,T1 = D(1,m1,T1)
                    L1.append(k1)
                    k1,m1,T1 = D01(m1,T1)
                    L1.append(k1)
                    i = i+2
                else:
                    k1,m1,T1 = D(C[i],m1,T1)
                    L1.append(k1)
                    i = i+1
            if L == L1:
                i = 0
                while i < len(Cr):
                    if i < len(Cr)-1 and Cr[i+1] == 0:
                        m1,T1 = S01(Lr[i]+1,m1,T1)
                        for j in range(Lr[i+1]):
                            m1,T1 = S(1,m1,T1)
                        i = i+2
                    else:
                        for j in range(Lr[i]+1):
                            m1,T1 = S(Cr[i],m1,T1)
                        i = i+1
                Pi.append((m1,T1))
        return Pi
    if s < 0 and s >= -(b-1)/2:
        m,T = rep(E)
        C = []
        for j in range(-int(2*s)):
            for i in range(a):
                C.append(B+s-i+j)
        L = []
        i = 0
        while i < len(C):
            if C[i] == 0:
                k,m,T = D(-1,m,T)
                L.append(k)
                k,m,T = D00(m,T)
                L.append(k)
                i = i+2
            else:
                k,m,T = D(C[i],m,T)
                L.append(k)
                i = i+1
        Cr = list(reversed(C)); Lr = list(reversed(L))                
        Pi1 = decomp((a,int(b+int(2*s))),E)
        for E1 in Pi1:
            m1,T1 = rep(E1)
            L1 = []
            i = 0
            while i < len(C):
                if C[i] == 0:
                    k1,m1,T1 = D(-1,m1,T1)
                    L1.append(k1)
                    k1,m1,T1 = D00(m1,T1)
                    L1.append(k1)
                    i = i+2
                else:
                    k1,m1,T1 = D(C[i],m1,T1)
                    L1.append(k1)
                    i = i+1
            if L == L1:
                i = 0
                while i < len(Cr):
                    if i < len(Cr)-1 and Cr[i+1] == 0:
                        for j in range(Lr[i]+1):
                            m1,T1 = S00(m1,T1)
                        for j in range(Lr[i+1]):
                            m1,T1 = S(-1,m1,T1)
                        i = i+2
                    else:
                        for j in range(Lr[i]+1):
                            m1,T1 = S(Cr[i],m1,T1)
                        i = i+1
                Pi.append((m1,T1))
        return Pi


###################################################################
# Speh times Arthur: irreducibility and first reducible points.
###################################################################

def Is_irred(s,ab,E):
    """
    Irreduciblity of Speh(a,b)||**{s} \rtimes rep(E).
    Input: (s,a,b,E)
        s (real number)
        ab (= (a,b), two positive integers)
        E (extended multi-segment)
    Output: True or False
        True <=> Speh(a,b)||**{s} \rtimes rep(E) is irreducible.
    """
    a = ab[0]; b = ab[1]
    Pi1 = soc(s,(a,b),E); Pi2 = soc(-s,(a,b),E)
    return len(Pi1) == 1 and len(Pi2) == 1 and Pi1 == Pi2

def FRP(ab,E):
    """
    First reducibility point for Speh(a,b)||**{s} \rtimes rep(E).
    Input: (s,a,b,E)
        ab (= (a,b), two positive integers)
        E (extended multi-segment)
    Output: s
        s (non-negative real number)
            the minimal non-negative real number  
            such that Speh(a,b)||**{s} \rtimes rep(E) is reducible.
    """
    a = ab[0]; b = ab[1]
    if (a-b - 2*E[0][0][0]) % 2 == 0:
        s = 0
    else: 
        s = 1/2
    while Is_irred(s,(a,b),E):
        s = s+1
    return s


###########################################################
# Bug report.
###########################################################

