#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
//This fuction should measure <Nhc Si^z Sj^z> <Nhc Si^+ Sj^-> <Nhc Si^- Sj^+>
// Where i<c<j
// Nh is a hole number operator
//------------------------------------------------------------------
void chargespinspinCorrelator(MPS psi, const SiteSet sites, int c){
	auto N = length(psi);
	printfln("i j c <Nhc Si^z Sj^z> <Nhc Si^+ Sj^-> <Nhc Si^- Sj^+>");
        for (int i=1;i<=N;i++)
	{ 
         for (int j=i;j<=N;j++)
	    {
	     if (i<c)
	     {   psi.position(i);
		auto Cz = psi(i)*sites.op("Sz",i);
		auto Cpm = psi(i)*sites.op("S+",i);
	        auto Cmp = psi(i)*sites.op("S-",i); 
		if ((j>i) && (j>c))
		 {
		   auto ir = commonIndex(psi(i),psi(i+1),"Link");
                   Cz *= dag(prime(prime(psi(i),"Site"),ir));
		   Cpm *= dag(prime(prime(psi(i),"Site"),ir));
		   Cmp *= dag(prime(prime(psi(i),"Site"),ir));
                  for (int k=i+1;k<c;k++){
			    Cz *= psi(k);
	 		    Cz *= dag(prime(psi(k),"Link"));
			    Cpm *= psi(k);
			    Cpm *= dag(prime(psi(k),"Link"));
			    Cmp *= psi(k);
			    Cmp *= dag(prime(psi(k),"Link"));
			    }	    
		Cz *= psi(c)*sites.op("Nh",c);
	        Cpm *= psi(c)*sites.op("Nh",c);
          	Cmp *= psi(c)*sites.op("Nh",c); 
	        Cz *= dag(prime(prime(psi(c),"Site"),"Link"));
	        Cpm *= dag(prime(prime(psi(c),"Site"),"Link"));
	        Cmp *= dag(prime(prime(psi(c),"Site"),"Link"));
         	for (int h=c+1;h<j;h++){
		     Cz *= psi(h);
		     Cz *= dag(prime(psi(h),"Link"));
		     Cpm *= psi(h);
		     Cpm *= dag(prime(psi(h),"Link"));
		     Cmp *= psi(h);
		     Cmp *= dag(prime(psi(h),"Link"));
		     }
		Cz *= psi(j)*sites.op("Sz",j);
		Cpm *= psi(j)*sites.op("S-",j);
	        Cmp *= psi(j)*sites.op("S+",j);
	        auto il = commonIndex(psi(j),psi(j-1),"Link");
		Cz *= dag(prime(prime(psi(j),"Site"),il));
		Cpm *= dag(prime(prime(psi(j),"Site"),il));
	        Cmp *= dag(prime(prime(psi(j),"Site"),il));
		} //closes if j>i and j>c
	printfln("",i," ",j," ",elt(Cz)," ",eltC(Cpm)," ",eltC(Cmp));	
	     } // closes if i<c
	         } // closes for (int j=i;j<=N;j++)
	         } // closes for (int i=1;i<=N;i++)
	         } // closes the function

