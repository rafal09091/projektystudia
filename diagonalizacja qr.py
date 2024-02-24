import sys
import numpy as np



dostabilnosci=0
kiedyzero=0.01
domyslnatolerancja=0.0001
maxliczbapowt=1000


def Haus(xx, n):
    x=xx.copy()
    y=x[n:,n:]
    Hs=np.eye(y.shape[0])
    v=y[:,0]
    c=dostabilnosci+(v.reshape(1,-1)@v.reshape(-1,1))
    v=v/np.sqrt(c)
    e=np.zeros(v.shape[1])
    e[0]=1.0
    w=v+e
    d=dostabilnosci+(w.reshape(1,-1)@w.reshape(-1,1))
    w=w/np.sqrt(d)
    Hs=Hs-2*w.reshape(-1,1)@w.reshape(1,-1)
    H=np.eye(x.shape[0])
    H[n:,n:]=Hs
    return H
    
def twoxtwo(x):
    D=x[0,0]*x[1,1]-x[0,1]*x[1,0]
    T=x[0,0]+x[1,1]
    delta=T*T-4*D
    lambda1=0.5*complex(T,np.sqrt(np.abs(delta)))
    lambda2=0.5*complex(T,-np.sqrt(np.abs(delta)))
    return lambda1,lambda2

def QR(x):
    Q=np.eye(x.shape[0])
    R=x.copy()    
    for i in range(x.shape[0]-1):
        Ht=Haus(R,i)
        Q=Ht@Q
        R=Ht@R
    return Q.T,R

def diag(x,tre = domyslnatolerancja, powt=maxliczbapowt):
    Q,R=QR(x)
    A=x.copy()
    e=np.max(np.abs(Q)-np.eye(Q.shape[0]))
    vectors=np.eye(x.shape[0])
    vectors=vectors
    licz=0
    while((e>tre) & (licz<powt)):
        Q,R=QR(A)
        A=R@Q
        vectors=vectors@Q
        e=np.max(np.abs(Q)-np.eye(Q.shape[0]))
        licz=licz+1
    return A,vectors.T

def eigenvectorS(x, n):
    y=x.copy()
    y=y[:(n+1),:(n+1)]
    v=np.zeros(x.shape[0])
    v[n]=1.0
    i=n-1
    for k in range(y.shape[0]):
        y[k,k]=y[k,k]-y[n,n]
    while(i>=0):
        assert (y[i,i]!=0), "nie da sie"
        s=0.0
        for k in range(i+1,n+1):
            s=s+v[k]*y[i,k]
        v[i]=-1*s/y[i,i]
        i=i-1
    return x[n,n],v

def eigenvectorC1(x, n):
    y=x.astype(complex)
    y=y[:(n+2),:(n+2)]
    #re=x[n,n]
    #im=x[n+1,n]
    tau,_=twoxtwo(y[n:n+2,n:n+2])
    v=np.zeros(x.shape[0]).astype(complex)
    v[n+1]=1.0
    v[n]=-x[n,n+1]/(x[n,n]-tau)
    i=n-1
    for k in range(y.shape[0]):
        y[k,k]=y[k,k]-tau
    while(i>=0):
        assert (y[i,i]!=0), "nie da sie"
        s=0.0
        for k in range(i+1,n+2):
            s=s+v[k]*y[i,k]
        v[i]=-1*s/y[i,i]
        i=i-1
    return tau,v


def eigenvectorC2(x, n):
    y=x.astype(complex)
    y=y[:(n+1),:(n+1)]
    #re=x[n,n]
    #im=x[n-1,n]
    _,tau=twoxtwo(y[n-1:n+1,n-1:n+1])
    v=np.zeros(x.shape[0]).astype(complex)
    v[n]=1.0
    v[n-1]=-x[n-1,n]/(x[n-1,n-1]-tau)
    i=n-2
    for k in range(y.shape[0]):
        y[k,k]=y[k,k]-tau
    while(i>=0):
        assert (y[i,i]!=0), "nie da sie"
        s=0.0
        for k in range(i+1,n+1):
            s=s+v[k]*y[i,k]
        v[i]=-1*s/y[i,i]
        i=i-1
    return tau,v


#funkcja do obliczania wektorow wlasnych numer n z macierzy x
def eigenvector(x,n):
    #gdy n jest rowne zero
    if (n==0):
        #czy chodzi o wartosc wlasna rzeczywista czy zespolona - pytam czy jeden element pod diagonala jest rowny zero
        if(np.abs(x[1,0])<kiedyzero):
            return eigenvectorS(x,n)
        else:
            return eigenvectorC1(x,n)
    #gdy n jest ostatnim wektorem
    elif (n==x.shape[0]-1):
        ##czy chodzi o wartosc wlasna rzeczywista czy zespolona - pytam czy jeden na lewo jest rowny zero
        if(np.abs(x[n,n-1])<kiedyzero):
            return eigenvectorS(x,n)
        else:
            return eigenvectorC2(x,n)
 #gdy n jest pomiedzy
    else:
        #czy chodzi o wartosc wlasna rzeczywista czy zespolona - pytam czy jeden element pod diagonala jest rowny zero lub czy jeden element na lewo jest rowny zero
        if((np.abs(x[n,n-1])<kiedyzero)&(np.abs(x[n+1,n])<kiedyzero)):
            return eigenvectorS(x,n)
        elif(np.abs(x[n,n-1])<kiedyzero):
            return eigenvectorC1(x,n)
        elif(np.abs(x[n+1,n])<kiedyzero):
            return eigenvectorC2(x,n)
        else:
            assert False, "nie da sie"
        


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

    print(f"wymiar (n): {n}")
    print("Macierz:")
    print(matrix)
    #domyslna tolerancja to dopuszalny blad, kryterium to bliskosc obrotu do bycia identycznoscia, maks liczba powtorzen
    D,VS=diag(matrix,domyslnatolerancja,maxliczbapowt)
    print(D)
    for i in range(D.shape[0]):
        wartosc,evp=eigenvector(D, i)
        ev= (VS.T)@evp
        print(f'wartosc wlasna {np.array([wartosc])} jej wektor wlasny {ev}' )



if __name__ == "__main__":
    main()
