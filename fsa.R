
# R Functions for Basic Frequency Stability Analysis
# W.J. Riley
# Hamilton Technical Services, Beaufort, SC 29907 USA
# Version 0.1.0
# May 16, 2020

# Packages required
library(allanvar) # For avar() and avari()
library(RobPer) # For TK95()
library(zoo) # For rollapply()

# Function to average phase data
pavg<-function(x,af=2)
{
  seed<-c(TRUE,rep(FALSE,af-1))
  cont<-rep(seed,ceiling(length(x)/af))[1:length(x)]
  return(x[which(cont)])
}

# Function to average frequency data
favg<-function(data,af=2)
{
  return(rollapply(data,width=af,by=af,FUN=mean))
}

# Function for phase to frequency conversion
ptof<-function(x,tau=1)
{
  return(diff(x)/tau)
}

# Function for frequency to phase conversion
ftop<-function(y,tau=1)
{
  return(diffinv(y)*tau)
}

# Function to generate power law noise
noise<-function(num,alpha,sigma,type=0,tau=1)
{
  z<-TK95(num,alpha)
  if(type==0) d<-padev(z,tau)
  else d<-fadev(z)
  z=(z/d)*sigma
  m=mean(z)
  z=z-m
  return(z)
}

# Function to find the noise type using the lag 1 ACF method
nid<-function(z)
{
  nD=0 # Difference order
  # Save original data
  zz<-z
  # Calc lag 1 autocorrelation r1
  r1=acf(z,1, "cor",F)
  r1=r1$acf[2]
  # Find d = r1/(1+r1)
  d=r1/(1+r1)
  # If d<0.25, must apply increment operator
  if(d>0.25)
  {
    while(d>=0.25)
    {
      # Take 1st differences
      z<-diff(z)
      nD=nD+1
      # Calc lag 1 autocorrelation r1
      r1=acf(z,1, "cor",F)
      r1=r1$acf[2]
      # Find d = r1/(1+r1)
      d=r1/(1+r1)
    }
  }
  # Calc alpha
  alpha=-2*d -2*nD +2
  # Restore original data
  z<-zz
  return (alpha)
}

# Function to show basic statistics for phase or frequency data
bs <- function(z,type=0,tau=1)
{
  print("Basic Statistics:", quote=FALSE)
  txt=paste("File =", deparse(substitute(z)))
  print(txt, quote=FALSE)
  if(type==0)
  {
    print("Type = Phase", quote=FALSE)
  }
  else
  {
    print("Type = Frequency", quote=FALSE)
  }
  txt=paste("Tau =", tau)
  print(txt, quote=FALSE)
  txt=paste("# Points =", length(z))
  print(txt, quote=FALSE)
  txt=paste("Max =", max(z))
  print(txt, quote=FALSE)
  txt=paste("Min =", min(z))
  print(txt, quote=FALSE)
  txt=paste("Span=", max(z)-min(z))
  print(txt, quote=FALSE)
  txt=paste("Mean =", mean(z))
  print(txt, quote=FALSE)
  txt=paste("Median =", median(z))
  print(txt, quote=FALSE)
  txt=paste("MAD =", mad(z))
  print(txt, quote=FALSE)
  txt=paste("Std Dev =", sqrt(var(z)))
  print(txt, quote=FALSE)
  if(type==0) # Phase data
  {
    txt=paste("Sigma =", padev(z,tau))
    print(txt, quote=FALSE)
    # txt=paste("nid =", nid(z))
    # print(txt, quote=FALSE)
    alpha=nid(z)
  }
  else # Freq data
  {
    txt=paste("Sigma =", fadev(z))
    print(txt, quote=FALSE)
    # txt=paste("nid =", nid(z))
    # print(txt, quote=FALSE)
    alpha=nid(z)-2
  }
  txt=paste("Alpha =", alpha )
  print(txt, quote=FALSE)
  if(alpha>1.5)
  {
    txt=paste("Noise = W PM")
    print(txt, quote=FALSE)
  }
  else if(alpha>0.5)
  {
    txt=paste("Noise = F PM")
    print(txt, quote=FALSE)
  }
  else if(alpha>-0.5)
  {
    txt=paste("Noise = W FM")
    print(txt, quote=FALSE)
  }
  else if(alpha>-1.5)
  {
    txt=paste("Noise = F FM")
    print(txt, quote=FALSE)
  }
  else
  {
    txt=paste("Noise = RW FM")
    print(txt, quote=FALSE)
  }
  plot(z)
}

# Function to count outliers in phase or frequency data
co <- function(z, limit=5)
{
  # Find MAD
  m=mad(z)
  # Count outliers
  n=sum(z<(-m*limit))+sum(z>(m*limit))
  return (n)
}

# Function to calculate the ADEV for phase data
padev <- function (x, tau=1)
{
  N=length(x)
  s=0
  for (i in 1:(N-2))
  {
    s = s + (x[i+2]-(2*x[i+1])+x[i])^2
  }
  av = s/(2*(tau^2)*(N-2))
  return (sqrt(av))
}

# Function to calculate the ADEV for frequency data
fadev <- function(y)
{
  N=length(y)
  s=0
  for (i in 1:(N-1))
  {
    s = s + (y[i+1]-y[i])^2
  }
  av=s/(2*(N-1))
  return (sqrt(av))
}

# Function to calculate the overlapping Allan deviation from phase data
poadev <- function(x, tau=1, m=1)
{
  N=length(x)
  s=0
  for(i in 1:(N-2*m))
  {
    s = s + (x[i+2*m]-2*x[i+m]+x[i])^2
  }
  s = s/(2*m^2*(N-2*m)*tau^2)
  return (sqrt(s))
}

# Function to calculate the overlapping Allan deviation from frequency data
foadev <- function(y, tau=1, af=1)
{
  x=ftop(y,tau)
  ad=poadev(x,tau,af)
  return (ad)
}

# Function to calculate the overlapping Allan deviation from phase or data
# over a range of octave averaging factors
adevrun <- function(z, type=0, tau=1)
{
  # If frequency data, convert it to phase data
  if(type==1)
  {
    x<-ftop(z)
  }
  else
  {
    x<-z
  }
  # Initializations
  N=length(x)
  af=1
  # Loop thru AFs up to limit, calculating ADEV
  # The maximum AF is floor(N/4)
  while(af<=floor(N/4))
  {
    ad=poadev(x,tau,af)
    print(paste0("AF= ",af,"  ADEV=",ad))
    af=af*2
  }
}

# Function to calculate Modified Allan deviation
# MVAR for phase data
# Argument tau is basic data sampling interval
# Each analysis tau is tau*m
# where argument m is averaging factor 1 to N/3
pmdev<-function(x,tau=1,m=1)
{
  N=length(x)
  mvar=0
  # Outer loop
  for(j in 1:(N-3*m+1))
  {
    s=0
    # Inner loop
    for(i in j:(j+m-1))
    {
      s=s+(x[i+(2*m)]-2*x[i+m]+x[i])
    }
    mvar=mvar+s^2
  }
  # Scaling
  mvar=mvar/(2*m^2*m^2*tau^2*(N-3*m+1))
  return (sqrt(mvar))
}

# Function to calculate the Hadamard deviation for phase data
phdev <- function (x, tau=1)
{
  N=length(x)
  s=0
  for (i in 1:(N-3))
  {
    s = s +(x[i+3] -3*x[i+2] +3*x[i+1] -x[i])^2
  }
  hv = s/(6*(tau^2)*(N-3))
  return (sqrt(hv))
}

# Function to calculate the Hadamard deviation for frequency data
fhdev <- function(y)
{
  N=length(y)
  s=0
  for (i in 1:(N-2))
  {
    s = s + (y[i+2] -2*y[i+1] +y[i])^2
  }
  hv=s/(6*(N-2))
  return (sqrt(hv))
}

# Find Theo1 per Howe and Peppler (2003)
# x = phase data vector (1 to N)
# tau = data sampling interval
# m = averaging factor (2 to N-1)
# m must be even
# Analysis tau = m*tau
# Stride = 0.75*m*tau
theo1<-function(x, tau=1, m=2)
{
  # Initializations
  N-length(x)
  t1=0

  # Outer sum
  for( i in 1:(N-m))
  {
    sum=0
    # Inner sum
    for( d in 0:((m/2)-1))
    {
      s=(1/((m/2)-d))*((x[i]-x[i-d+(m/2)]+x[i+m]-x[i+d+(m/2)])^2)
      sum=sum+s
    }
    t1=t1+sum
  }

  # Scaling factor
  t1=t1/(0.75*(N-m)*(m*tau)^2)

  # Return Theo1 deviation
  return (sqrt(t1))
}

# Function to calculate and plot a power spectral density
psd <- function(z, span=10, logx=TRUE, logy=TRUE, title="PSD Plot")
{
  s<-spectrum(z,span)
  freq<-s$freq
  psd<-2*s$spec

  if( logx==FALSE & logy==FALSE)
  {
    plot(freq,psd,type="l",main=title)
  }
  else if(logx==FALSE & logy==TRUE)
  {
    plot(freq,log10(psd),type="l",main=title)
  }
  else if(logx==TRUE & logy==FALSE)
  {
    plot(log10(freq),psd,type="l",main=title)
  }
  else(logx==TRUE & logy==TRUE)
  {
    plot(log10(freq),log10(psd),type="l",main=title)
  }
}

