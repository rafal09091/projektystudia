
import copy
import sys

def der(B):
    return [i*j for i,j in enumerate(B)][1:]

def poly(w,x):
    y=1.0
    s=0.0
    for a in w:
        s=s+y*a
        y=y*x
    return s

def newton(w, x0=complex(0.0001,0.0001), tr=0.0001, bound=1000000):
    y=x0
    wp=der(w)
    e=10000.0
    while(e>tr):
        z=y-poly(w,y)/poly(wp,y)
        assert (abs(z)<bound), "nie ma zbieznosci"
        e=abs(z-y)
        y=z
    return y


def horner(w, x0):
    r = copy.deepcopy(w)
    v = copy.deepcopy(w)
    l = len(r) - 1
    while (l > 0):
        v[l - 1] = r[l]
        r[l - 1] = r[l - 1] + x0 * v[l - 1]
        r[l] = 0
        l = l - 1
    return v[:-1]

def findall(w,x0=complex(0.01,-0.01), tr=0.0001, bound=1000000):
    zeros=[]
    v=copy.deepcopy(w)
    prev=x0
    while (len(v)>2):
        y=newton(v,prev,tr, bound)
        zeros.append(y)
        v=horner(v,y)
        prev=x0
    zeros.append(-v[0]/v[1])
    return zeros
# findall podajemy wielomian [a_0, ... , a_n], punkt startowy, warunek zatrzymania
#print(findall([1.0,0.0,1.0]))
#print(findall([-1.0,0.0,1.0]))
#print(findall([0.0,-1.0,0.0,1.0],complex(0.5,0.5)))
name=sys.argv[1]
string=''

with open (name, "r") as f:
    lines = [
        line.strip("\n").split(",")
        for line in f.readlines()
    ]
    lines=lines[1:]
polynomialfromfile=list(map(lambda x: complex(float(x[0]),float(x[1])) if len(x)==2 else complex(float(x[0]),0.0) ,lines))
for i,j in enumerate(polynomialfromfile):
    string+= str(j)+' x^'+str(i)
    if (i!=len(polynomialfromfile)-1):
        string+=' + '
abss= list(map(lambda x: abs(x), polynomialfromfile))
ogr= len(abss)*max(1,sum(abss[0:-1])/abss[-1])
print('ograniczenie to: ' )
print(ogr)
print('wczytany wielomian to: '+string)
print('zera wielomianu to: ')
print(findall(polynomialfromfile,bound=ogr*1000))
