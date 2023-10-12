n=200
p=(1:99)/100
p
mcchi=200000
chis=c()
for (i in 1:mcchi)
{
U=rnorm(4,0,1)
chis[i]=sum(U^2)
}
pt=c() ## p wartosci
X=c()  
mc=10000
Q=c()
moce=c()
for (k in 1:99)
{
for (j in 1:mc)
{
  N=c(0,0,0,0,0,0)
  for (i in 1:n)
  {
    D=rbinom(1,1,0.5)
    if (D<0.5)
    {
      X[i]=rbinom(1,5,0.5)
      N[X[i]+1]=N[X[i]+1]+1
    }
    else
    {
      X[i]=rbinom(1,5,p[k])
      N[X[i]+1]=N[X[i]+1]+1
    }
  }
  phat=sum(X)/5/n
  g=c()
  for (i in (1:6))
  {
    g[i]=factorial(5)/factorial(i-1)/factorial(5-i+1)*phat^(i-1)*(1-phat)^(5-i+1)
  }
  Q[j]=sum((N-n*g)^2/(n*g))
  pt[j]=sum(chis>Q[j])/mcchi
}
moce[k]=mean(pt<0.05)
print(k)
}
plot(moce~p)
abline(h=0.05)
