#include "itensor/all.h"
#include "itensor/util/print_macro.h"

//Nh is a hole number operator that has to be added manualy

using namespace itensor;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////computing charge spin spin correlation
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void chargespinspinCorrelator(MPS psi, const SiteSet sites, int c){
	auto N = length(psi);
	printfln("i j <Nhc Si^z Sj^z> <Nhc Si^+ Sj^-> <Nhc Si^- Sj^+>");
        for (int i = 1;i<=N;i++)
	{if(c>i)
	{for(int j=i;j<c;j++)
                {  //----------------------------i=j c----------------------------------
		if(j==i){psi.position(i);
                auto Cz = psi(i)*sites.op("Sz",i)*prime(sites.op("Sz",i));
                auto Cpm = psi(i)*sites.op("S+",i)*prime(sites.op("S-",i));
                auto Cmp = psi(i)*sites.op("S-",i)*prime(sites.op("S+",i));
                auto ir = commonIndex(psi(i),psi(i+1),"Link");
		   Cz *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ir));
                   Cpm *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ir));
                   Cmp *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ir));
		   for (int h=j+1;h<c;h++){
                     Cz *= psi(h);
                     Cz *= dag(prime(psi(h),"Link"));
                     Cpm *= psi(h);
                     Cpm *= dag(prime(psi(h),"Link"));
                     Cmp *= psi(h);
                     Cmp *= dag(prime(psi(h),"Link"));
                     }
                Cz *= psi(c)*sites.op("Nh",c);
                Cpm *= psi(c)*sites.op("Nh",c);
                Cmp *= psi(c)*sites.op("Nh",c);
                auto ill = commonIndex(psi(c),psi(c-1),"Link");
                Cz *= dag(prime(prime(psi(c),"Site"),ill));
                Cpm *= dag(prime(prime(psi(c),"Site"),ill));
                Cmp *= dag(prime(prime(psi(c),"Site"),ill));
		printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp)); 
		}//closes if i==j
	 else{//-------------------------------------------------i j c-------------------------
		 psi.position(i);
                auto Cz = psi(i)*sites.op("Sz",i);
                auto Cpm = psi(i)*sites.op("S+",i);
                auto Cmp = psi(i)*sites.op("S-",i);
               auto ir = commonIndex(psi(i),psi(i+1),"Link");

                   Cz *= dag(prime(prime(psi(i),"Site"),ir));
                   Cpm *= dag(prime(prime(psi(i),"Site"),ir));
                   Cmp *= dag(prime(prime(psi(i),"Site"),ir));
	   for (int k=i+1;k<j;k++){
                            Cz *= psi(k);
                            Cz *= dag(prime(psi(k),"Link"));
                            Cpm *= psi(k);
                            Cpm *= dag(prime(psi(k),"Link"));
                            Cmp *= psi(k);
                            Cmp *= dag(prime(psi(k),"Link"));
                            }
	   Cz *= psi(j)*sites.op("Sz",j);
                Cpm *= psi(j)*sites.op("S-",j);
                Cmp *= psi(j)*sites.op("S+",j);
                Cz *= dag(prime(prime(psi(j),"Site"),"Link"));
                Cpm *= dag(prime(prime(psi(j),"Site"),"Link"));
                Cmp *= dag(prime(prime(psi(j),"Site"),"Link"));
		for (int h=j+1;h<c;h++){
                     Cz *= psi(h);
                     Cz *= dag(prime(psi(h),"Link"));
                     Cpm *= psi(h);
                     Cpm *= dag(prime(psi(h),"Link"));
                     Cmp *= psi(h);
                     Cmp *= dag(prime(psi(h),"Link"));
                     }
		Cz *= psi(c)*sites.op("Nh",c);
                Cpm *= psi(c)*sites.op("Nh",c);
                Cmp *= psi(c)*sites.op("Nh",c);
                auto ill = commonIndex(psi(c),psi(c-1),"Link");
                Cz *= dag(prime(prime(psi(c),"Site"),ill));
                Cpm *= dag(prime(prime(psi(c),"Site"),ill));
                Cmp *= dag(prime(prime(psi(c),"Site"),ill));
printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
	 }//closes else

}// closes for j=i to c
	 //------------------------------------------i c j-------------------------------------------------
		for (int j=c+1;j<=N;j++)
	{psi.position(i);
                auto Cz = psi(i)*sites.op("Sz",i);
                auto Cpm = psi(i)*sites.op("S+",i);
                auto Cmp = psi(i)*sites.op("S-",i);
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
	        auto ill = commonIndex(psi(j),psi(j-1),"Link");
		Cz *= dag(prime(prime(psi(j),"Site"),ill));
		Cpm *= dag(prime(prime(psi(j),"Site"),ill));
	        Cmp *= dag(prime(prime(psi(j),"Site"),ill));
	printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));	
	} //closes for j=c+1
} //closes if i<c
else if (i>c)
{ //-----------------------------------------------c i=j-------------------------------------
	for(int j=i;j<=N;j++)
                {
                if(j==i){psi.position(c);
                auto Cz = psi(c)*sites.op("Nh",c);
                auto Cpm = psi(c)*sites.op("Nh",c);
                auto Cmp = psi(c)*sites.op("Nh",c);
                auto ir = commonIndex(psi(c),psi(c+1),"Link");
                   Cz *= dag(prime(prime(psi(c),"Site"),ir));
                   Cpm *= dag(prime(prime(psi(c),"Site"),ir));
                   Cmp *= dag(prime(prime(psi(c),"Site"),ir));
                   for (int h=c+1;h<i;h++){
                     Cz *= psi(h);
                     Cz *= dag(prime(psi(h),"Link"));
                     Cpm *= psi(h);
                     Cpm *= dag(prime(psi(h),"Link"));
                     Cmp *= psi(h);
                     Cmp *= dag(prime(psi(h),"Link"));
                     }
                Cz *= psi(i)*sites.op("Sz",i)*prime(sites.op("Sz",i));
                Cpm *= psi(i)*sites.op("S+",i)*prime(sites.op("S-",i));
                Cmp *= psi(i)*sites.op("S-",i)*prime(sites.op("S+",i));
                auto ill = commonIndex(psi(i),psi(i-1),"Link");
                Cz *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ill));
                Cpm *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ill));
                Cmp *= dag(prime(prime(prime(psi(i),"Site"),"Site"),ill));
                printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
                }//closes if i==j
		else
		{ //---------------------------------------------c i j--------------------------------------------------
			psi.position(c);
                auto Cz = psi(c)*sites.op("Nh",c);
                auto Cpm = psi(c)*sites.op("Nh",c);
                auto Cmp = psi(c)*sites.op("Nh",c);
               auto ir = commonIndex(psi(c),psi(c+1),"Link");

                   Cz *= dag(prime(prime(psi(c),"Site"),ir));
                   Cpm *= dag(prime(prime(psi(c),"Site"),ir));
                   Cmp *= dag(prime(prime(psi(c),"Site"),ir));
           for (int k=c+1;k<i;k++){
                            Cz *= psi(k);
                            Cz *= dag(prime(psi(k),"Link"));
                            Cpm *= psi(k);
                            Cpm *= dag(prime(psi(k),"Link"));
                            Cmp *= psi(k);
                            Cmp *= dag(prime(psi(k),"Link"));
                            }
           Cz *= psi(i)*sites.op("Sz",i);
                Cpm *= psi(i)*sites.op("S+",i);
                Cmp *= psi(i)*sites.op("S-",i);
                Cz *= dag(prime(prime(psi(i),"Site"),"Link"));
                Cpm *= dag(prime(prime(psi(i),"Site"),"Link"));
                Cmp *= dag(prime(prime(psi(i),"Site"),"Link"));
                for (int h=i+1;h<j;h++){
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
                auto ill = commonIndex(psi(j),psi(j-1),"Link");
                Cz *= dag(prime(prime(psi(j),"Site"),ill));
                Cpm *= dag(prime(prime(psi(j),"Site"),ill));
                Cmp *= dag(prime(prime(psi(j),"Site"),ill));
printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
		} //closes else
		} //closes for j

} // closes if i>c
} // closes for i 
} // closes the function

