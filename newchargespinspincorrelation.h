#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;


void new_chragespinCorrelator(MPS psi, const SiteSet sites, int c){
auto N = length(psi);
printfln("i j <Nhc Si^z Sj^z> <Nhc Si^+ Sj^-> <Nhc Si^- Sj^+>");
for (int i=1;i<=N;i++)
{
  for(int j=1;j<=i;j++)
  {
    auto Cz = psi(1);
    auto Cpm = psi(1);
    auto Cmp = psi(1);
if ((i!=j) && (j!=c) && (i!=c))
   {
     if (j==1) {Cz *= sites.op("Sz",1)*dag(prime(psi(1),"Site"));
                Cpm *= sites.op("S+",1)*dag(prime(psi(1),"Site"));
                Cmp *= sites.op("S-",1)*dag(prime(psi(1),"Site"));}
     else {Cz *= dag(psi(1));
           Cpm *= dag(psi(1));
           Cmp *= dag(psi(1));}
     for (int n=2;n<=N;n++)
     {
       Cz *= psi(n);
       Cpm *= psi(n);
       Cmp *= psi(n);
       if (j==n) { Cz *= sites.op("Sz",n)*dag(prime(psi(n),"Site"));
	           Cpm *= sites.op("S+",n)*dag(prime(psi(n),"Site"));
	           Cmp *= sites.op("S-",n)*dag(prime(psi(n),"Site"));}
       else if (i==n) {Cz *= sites.op("Sz",n)*dag(prime(psi(n),"Site"));
	              Cpm *= sites.op("S-",n)*dag(prime(psi(n),"Site"));
		      Cmp *= sites.op("S+",n)*dag(prime(psi(n),"Site"));}
       else if (c==n) {Cz *= sites.op("Nh",n)*dag(prime(psi(n),"Site")); 
         	       Cpm *= sites.op("Nh",n)*dag(prime(psi(n),"Site"));
		       Cmp *= sites.op("Nh",n)*dag(prime(psi(n),"Site"));}
       else   {Cz *= dag(psi(n));
	       Cpm *= dag(psi(n));
	       Cmp *= dag(psi(n));}
     } //closes for(n=2 to N)
     printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
   } //closes if ((i!=j) && (j!=c) && (i!=c))
else if ((i==j) && (i!=c))
{ if (j==1) {Cz *= sites.op("Sz",1)*prime(sites.op("Sz",1))*dag(prime(prime(psi(1),"Site"),"Site"));
       	     Cpm *= sites.op("S+",1)*prime(sites.op("S-",1))*dag(prime(prime(psi(1),"Site"),"Site"));
	     Cmp *= sites.op("S-",1)*prime(sites.op("S+",1))*dag(prime(prime(psi(1),"Site"),"Site"));}
else {Cz *= dag(psi(1));    
	Cpm *= dag(psi(1));
	Cmp *= dag(psi(1));}
for (int n=2;n<=N;n++)      
{
	Cz *= psi(n);
	Cpm *= psi(n);
	Cmp *= psi(n);
        if (i==n){ Cz *= sites.op("Sz",n)*prime(sites.op("Sz",n))*dag(prime(prime(psi(n),"Site"),"Site"));
		Cpm *= sites.op("S+",n)*prime(sites.op("S-",n))*dag(prime(prime(psi(n),"Site"),"Site"));
		Cmp *= sites.op("S-",n)*prime(sites.op("S+",n))*dag(prime(prime(psi(n),"Site"),"Site"));}
	else if (c==n) {Cz *= sites.op("Nh",n)*dag(prime(psi(n),"Site"));
	                Cpm *= sites.op("Nh",n)*dag(prime(psi(n),"Site"));
	                Cmp *= sites.op("Nh",n)*dag(prime(psi(n),"Site"));}
	else {Cz *= dag(psi(n));
	      Cpm *= dag(psi(n));
	      Cmp *= dag(psi(n));}
	} //closes for (n=2 to N)
printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
 } //closes else if((i==j) && (i!=c))
  } //closes for(j=1 to i)

} //closes for(i=1 to N)

} //closes the function
