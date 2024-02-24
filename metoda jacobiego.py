import sys
import numpy as np
import random
from time import time 

def timer_func(func): 
    def wrap_func(*args, **kwargs): 
        t1 = time() 
        result = func(*args, **kwargs) 
        t2 = time() 
        print(f'Czas: {func.__name__!r} executed in {(t2-t1):.4f}s') 
        return result 
    return wrap_func 

#mnozenie macierzy przez macierz obrotu jakobiego R
def multLO(R,X,p,q):
    Y=X.copy()
    apq= R[p,q]
    aqp= R[q,p]
    app= R[p,p]
    aqq= R[q,q]
    v1=Y[:,p].copy()
    v2=Y[:,q].copy()
    Y[:,p]=app*v1+v2*aqp
    Y[:,q]=apq*v1+v2*aqq
    return Y
def multRO(R,X,p,q):
    Y=X.copy()
    apq= R[p,q]
    aqp= R[q,p]
    app= R[p,p]
    aqq= R[q,q]
    v1=Y[p].copy()
    v2=Y[q].copy()
    Y[p]=app*v1+v2*apq
    Y[q]=aqp*v1+v2*aqq
    return Y



def maks (x):
    y=x.copy()
    for i in range(x.shape[0]):
        y[i,i]=0.0
    ind = np.unravel_index(np.argmax(np.abs(y), axis=None), y.shape)
    return ind, x[ind]
def createR (x):
    R=np.eye(x.shape[0])
    ind, apq = maks(x)
    p=ind[0]
    q=ind[1]
    aqp= x[q,p]
    app= x[p,p]
    aqq= x[q,q]
    if (app==aqq):
        theta=0.7853981633974483
    else:
        theta=np.arctan(2*apq/(app-aqq))/2
    c=np.cos(theta)
    s=np.sin(theta)
    R[p,p]=c
    R[q,q]=c
    R[p,q]=-s
    R[q,p]=s
    return R, np.abs(apq),p,q
def createRfromIndex (x,p,q):
    R=np.eye(x.shape[0])
    apq= x[p,q]
    aqp= x[q,p]
    app= x[p,p]
    aqq= x[q,q]
    if (app==aqq):
        theta=0.7853981633974483
    else:
        theta=np.arctan(2*apq/(app-aqq))/2
    c=np.cos(theta)
    s=np.sin(theta)
    R[p,p]=c
    R[q,q]=c
    R[p,q]=-s
    R[q,p]=s
    return apq,R
    
    
@timer_func    
def iter(x, tr=0.001, powt= 10000):
    _,e= maks(np.abs(x))
    n=0
    z=x.copy()
    #RK to skumulowana macierz obrotow
    RK=np.eye(x.shape[0])
    while((e>tr) & (n<powt)):
        #zmienilem to teraz createR zwraca mi tez ile wynosilo maksimum poza przekatna 
        # tworze sobie nowa macierz obrotu odpowiednio wykonuje mnozenie 
        R,e,p,q=createR(z)
        #RK=RK@R
        RK=multRO(R,RK,p,q)
        #z= (R.T)@z@R
        z= multLO(R.T,multRO(R,z,p,q),p,q)       
        n=n+1
        #print(z)
    print (f'wykonano {n} krokow')
    print(z)
    return z,RK
@timer_func    
def iterCyclic(x, tr=0.001, powt= 10000):
    z=x.copy()
    sall=np.sum(z*z)
    sdiag=np.trace(z*z)
    e=sall-sdiag
    n=0 
    RK=np.eye(x.shape[0])
    while((e>tr) & (n<powt)):
        eloc=0.0
        #iteracja po elementach poza diagonala
        for i in range(1,x.shape[0]):
            for j in range(0,i):
                aij,R=createRfromIndex(z,i,j)
                #RK=RK@R
                RK=multRO(R,RK,i,j)
                #z= (R.T)@z@R
                z= multLO(R.T,multRO(R,z,i,j),i,j) 
                n=n+1
                #print(z)
        sdiag=np.trace(z*z)       
        e=sall-sdiag
    print (f'wykonano {n} krokow')
    print(z)
    return z,RK
@timer_func
def iterRandom(x, tr=0.001, powt= 10000):
    n=0
    z=x.copy()
    sall=np.sum(z*z)
    sdiag=np.trace(z*z)
    e=sall-sdiag
    RK=np.eye(x.shape[0])
    while((e>tr) & (n<powt)):
        eloc=0.0
        i=random.randint(1,x.shape[0]-1)
        j=random.randint(0,i-1)
        aij,R=createRfromIndex(z,i,j)
        #RK=RK@R
        RK=multRO(R,RK,i,j)
        #z= (R.T)@z@R
        z= multLO(R.T,multRO(R,z,i,j),i,j) 
        n=n+1
        #print(z)
        sdiag=np.trace(z*z)       
        e=sall-sdiag
    print (f'wykonano {n} krokow')
    print(z)
    return z,RK

    
    


def read_matrix_from_file(filename):
    with open(filename, 'r') as file:
        

        n = int(file.readline().strip())

        matrix = []
        for _ in range(n):
            row = list(map(float, file.readline().strip().split()))
            matrix.append(row)

    return n, np.array(matrix)

def main():
    np.set_printoptions(suppress=True)
    if len(sys.argv) != 2:
        print("Usage: python main.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    n, matrix = read_matrix_from_file(filename)

    print(f"Dimension (n): {n}")
    print("Macierz:")
    print(matrix)
    print("metoda max")
    # tr to akceptowalny blad powt to maksymalna liczba powtorzen
    S1,RK1=iter(matrix,tr=0.001, powt= 10000)
    print("   ")
    for i in range(matrix.shape[0]):
        print(f'wartosc wlasna {S1[i,i]}')
        print(f'wektor wlasny {RK1[:,i]}')
    print("metoda cykliczna")
    # tr to akceptowalny blad powt to maksymalna liczba powtorzen
    S2,RK2=iterCyclic(matrix,tr=0.001, powt= 10000)
    print("   ")
    for i in range(matrix.shape[0]): 
         print(f'wartosc wlasna {S2[i,i]}')
         print(f'wektor wlasny {RK2[:,i]}')
    print("metoda losowa")
    # tr to akceptowalny blad powt to maksymalna liczba powtorzen
    S3,RK3=iterRandom(matrix,tr=0.001, powt= 10000)
    print("   ")
    for i in range(matrix.shape[0]):
        print(f'wartosc wlasna {S3[i,i]}')
        print(f'wektor wlasny {RK3[:,i]}')

        

    



if __name__ == "__main__":
    main()
