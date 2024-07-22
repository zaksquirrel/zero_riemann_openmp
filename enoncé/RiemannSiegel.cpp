#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <cassert>

/*************************************************************************
* *


This code computes the number of zeros on the critical line of the Zeta function.
https://en.wikipedia.org/wiki/Riemann_zeta_function 

This is one of the most important and non resolved problem in mathematics : https://www.science-et-vie.com/article-magazine/voici-les-7-plus-grands-problemes-de-mathematiques-jamais-resolus

This problem has been the subject of one of the most important distributed computing project in the cloud : more than 10000 machines during 2 years. 
They used this algorithm: very well optimized.
This project failed, bitten by a smaller team that used a far better algorithm. 
The code is based on the Thesis of Glendon Ralph Pugh (1992) : https://web.viu.ca/pughg/thesis.d/masters.thesis.pdf

We can optimize the code in numerous ways, and obviously parallelize it. 

Remark: we do not compute the zeros: we count them to check that they are on the Riemann Line.
Remark: Andrew Odlyzko created a method that is far more efficient but too complex to be the subject of an algorithmetical tuning exercice. 

The exercise is to sample a region on the critical line to count how many times the function changes sign, so that there is at least 1 zero between 2 sampling points.
Here we use a constant sampling but you can recode entirely the way to proceed.

Only a correct (right) count matters, and the performance.

compile g++ RiemannSiegel.cpp -O -o RiemannSiegel
--------------------------------------------------------------------------
./RiemannSiegel 10 1000 100
I found 10142 Zeros in 3.459 seconds     # OK 
--------------------------------------------------------------------------
./RiemannSiegel 10 10000 10 
I found 10142 Zeros in 0.376 seconds     # OK
--------------------------------------------------------------------------
./RiemannSiegel 10 100000 10
I found 137931 Zeros in 6.934 seconds    # INCORRECT
--------------------------------------------------------------------------
./RiemannSiegel 10 100000 100
I found 138069 Zeros in 56.035 seconds   # OK
--------------------------------------------------------------------------
RiemannSiegel 10 1000000     need to find : 1747146     zeros
RiemannSiegel 10 10000000    need to find : 21136125    zeros
RiemannSiegel 10 100000000   need to find : 248888025   zeros
RiemannSiegel 10 1000000000  need to find : 2846548032  zeros
RiemannSiegel 10 10000000000 need to find : 32130158315 zeros


The more regions you validate and with the best timing, the more points you get.

The official world record of the zeros computed is 10^13 but with some FFTs and the method from Odlyzsko.
Compute time 1 year-core so an algortihm 10000*2*40 times more efficient than ZetaGrid's one. 

* *
*************************************************************************/

typedef unsigned long      ui32;
typedef unsigned long long ui64;

double dml_micros()
{
        static struct timezone tz;
        static struct timeval  tv;
        gettimeofday(&tv,&tz);
        return((tv.tv_sec*1000000.0)+tv.tv_usec);
}

int even(int n)
{
	if (n%2 == 0) return(1);
	else          return(-1);
}

double theta(double t)
{
	const double pi = 3.1415926535897932385;
	return(t/2.0*log(t/2.0/pi) - t/2.0 - pi/8.0 + 1.0/48.0/t + 7.0/5760.0/pow(t,3.0) + 31.0/80640.0/powl(t,5.0) +127.0/430080.0/powl(t,7.0)+511.0/1216512.0/powl(t,9.0));
	//https://oeis.org/A282898  // numerators
	//https://oeis.org/A114721  // denominators
}


double C(int n, double z){
	if (n==0)
		  return(.38268343236508977173 * pow(z, 0.0)
			+.43724046807752044936 * pow(z, 2.0)
			+.13237657548034352332 * pow(z, 4.0)
			-.01360502604767418865 * pow(z, 6.0)
			-.01356762197010358089 * pow(z, 8.0)
			-.00162372532314446528 * pow(z,10.0)
			+.00029705353733379691 * pow(z,12.0)
			+.00007943300879521470 * pow(z,14.0)
			+.00000046556124614505 * pow(z,16.0)
			-.00000143272516309551 * pow(z,18.0)
			-.00000010354847112313 * pow(z,20.0)
			+.00000001235792708386 * pow(z,22.0)
			+.00000000178810838580 * pow(z,24.0)
			-.00000000003391414390 * pow(z,26.0)
			-.00000000001632663390 * pow(z,28.0)
			-.00000000000037851093 * pow(z,30.0)
			+.00000000000009327423 * pow(z,32.0)
			+.00000000000000522184 * pow(z,34.0)
			-.00000000000000033507 * pow(z,36.0)
			-.00000000000000003412 * pow(z,38.0)
			+.00000000000000000058 * pow(z,40.0)
			+.00000000000000000015 * pow(z,42.0));
	else if (n==1)
		 return(-.02682510262837534703 * pow(z, 1.0)
			+.01378477342635185305 * pow(z, 3.0)
			+.03849125048223508223 * pow(z, 5.0)
			+.00987106629906207647 * pow(z, 7.0)
			-.00331075976085840433 * pow(z, 9.0)
			-.00146478085779541508 * pow(z,11.0)
			-.00001320794062487696 * pow(z,13.0)
			+.00005922748701847141 * pow(z,15.0)
			+.00000598024258537345 * pow(z,17.0)
			-.00000096413224561698 * pow(z,19.0)
			-.00000018334733722714 * pow(z,21.0)
			+.00000000446708756272 * pow(z,23.0)
			+.00000000270963508218 * pow(z,25.0)
			+.00000000007785288654 * pow(z,27.0)
			-.00000000002343762601 * pow(z,29.0)
			-.00000000000158301728 * pow(z,31.0)
			+.00000000000012119942 * pow(z,33.0)
			+.00000000000001458378 * pow(z,35.0)
			-.00000000000000028786 * pow(z,37.0)
			-.00000000000000008663 * pow(z,39.0)
			-.00000000000000000084 * pow(z,41.0)
			+.00000000000000000036 * pow(z,43.0)
			+.00000000000000000001 * pow(z,45.0));
else if (n==2)
		 return(+.00518854283029316849 * pow(z, 0.0)
			+.00030946583880634746 * pow(z, 2.0)
			-.01133594107822937338 * pow(z, 4.0)
			+.00223304574195814477 * pow(z, 6.0)
			+.00519663740886233021 * pow(z, 8.0)
			+.00034399144076208337 * pow(z,10.0)
			-.00059106484274705828 * pow(z,12.0)
			-.00010229972547935857 * pow(z,14.0)
			+.00002088839221699276 * pow(z,16.0)
			+.00000592766549309654 * pow(z,18.0)
			-.00000016423838362436 * pow(z,20.0)
			-.00000015161199700941 * pow(z,22.0)
			-.00000000590780369821 * pow(z,24.0)
			+.00000000209115148595 * pow(z,26.0)
			+.00000000017815649583 * pow(z,28.0)
			-.00000000001616407246 * pow(z,30.0)
			-.00000000000238069625 * pow(z,32.0)
			+.00000000000005398265 * pow(z,34.0)
			+.00000000000001975014 * pow(z,36.0)
			+.00000000000000023333 * pow(z,38.0)
			-.00000000000000011188 * pow(z,40.0)
			-.00000000000000000416 * pow(z,42.0)
			+.00000000000000000044 * pow(z,44.0)
			+.00000000000000000003 * pow(z,46.0));
else if (n==3)
		 return(-.00133971609071945690 * pow(z, 1.0)
			+.00374421513637939370 * pow(z, 3.0)
			-.00133031789193214681 * pow(z, 5.0)
			-.00226546607654717871 * pow(z, 7.0)
			+.00095484999985067304 * pow(z, 9.0)
			+.00060100384589636039 * pow(z,11.0)
			-.00010128858286776622 * pow(z,13.0)
			-.00006865733449299826 * pow(z,15.0)
			+.00000059853667915386 * pow(z,17.0)
			+.00000333165985123995 * pow(z,19.0)
			+.00000021919289102435 * pow(z,21.0)
			-.00000007890884245681 * pow(z,23.0)
			-.00000000941468508130 * pow(z,25.0)
			+.00000000095701162109 * pow(z,27.0)
			+.00000000018763137453 * pow(z,29.0)
			-.00000000000443783768 * pow(z,31.0)
			-.00000000000224267385 * pow(z,33.0)
			-.00000000000003627687 * pow(z,35.0)
			+.00000000000001763981 * pow(z,37.0)
			+.00000000000000079608 * pow(z,39.0)
			-.00000000000000009420 * pow(z,41.0)
			-.00000000000000000713 * pow(z,43.0)
			+.00000000000000000033 * pow(z,45.0)
			+.00000000000000000004 * pow(z,47.0));
else
		 return(+.00046483389361763382 * pow(z, 0.0)
			-.00100566073653404708 * pow(z, 2.0)
			+.00024044856573725793 * pow(z, 4.0)
			+.00102830861497023219 * pow(z, 6.0)
			-.00076578610717556442 * pow(z, 8.0)
			-.00020365286803084818 * pow(z,10.0)
			+.00023212290491068728 * pow(z,12.0)
			+.00003260214424386520 * pow(z,14.0)
			-.00002557906251794953 * pow(z,16.0)
			-.00000410746443891574 * pow(z,18.0)
			+.00000117811136403713 * pow(z,20.0)
			+.00000024456561422485 * pow(z,22.0)
			-.00000002391582476734 * pow(z,24.0)
			-.00000000750521420704 * pow(z,26.0)
			+.00000000013312279416 * pow(z,28.0)
			+.00000000013440626754 * pow(z,30.0)
			+.00000000000351377004 * pow(z,32.0)
			-.00000000000151915445 * pow(z,34.0)
			-.00000000000008915418 * pow(z,36.0)
			+.00000000000001119589 * pow(z,38.0)
			+.00000000000000105160 * pow(z,40.0)
			-.00000000000000005179 * pow(z,42.0)
			-.00000000000000000807 * pow(z,44.0)
			+.00000000000000000011 * pow(z,46.0)
			+.00000000000000000004 * pow(z,48.0));
}

double Z(double t, int n)
//*************************************************************************
// Riemann-Siegel Z(t) function implemented per the Riemenn Siegel formula.
// See http://mathworld.wolfram.com/Riemann-SiegelFormula.html for details
//*************************************************************************
{
	double p; /* fractional part of sqrt(t/(2.0*pi))*/
	double C(int,double); /* coefficient of (2*pi/t)^(k*0.5) */
	const double pi = 3.1415926535897932385; 
	int N = sqrt(t/(2.0 * pi)); 
	p = sqrt(t/(2.0 * pi)) - N; 
	double tt = theta(t); 
	double ZZ = 0.0; 
	for (int j=1;j <= N;j++) {
		ZZ = ZZ + 1.0/sqrt((double) j ) * cos(fmod(tt -t*log((double) j),2.0*pi));
	} 
	ZZ = 2.0 * ZZ; 
	double R  = 0.0; 
	for (int k=0;k <= n;k++) {
		R = R + C(k,2.0*p-1.0) * pow(2.0*pi/t, ((double) k)*0.5);
	} 
	R = even(N-1) * pow(2.0 * pi / t,0.25) * R; 
	return(ZZ + R);
}

/*
	Code to compute Zeta(t) with high precision
	Only works in IEEE 754 for t<1000
	This can help to validate that the Riemann Siegel function for small values but since we are mainly interrested to the behavior for large values of t,
	the best method is to compute zeros that are known
	As you may observe, the accuracy of Z(t) gets better with large values of t until being limited by the IEEE 754 norm and the double format. 
*/
std::complex <double> test_zerod(const double zero,const int N)
{
        std::complex <double> un(1.0,0);  
        std::complex <double> deux(2.0,0); 
        std::complex <double> c1(0.5,zero);
        std::complex <double> sum1(0.0,0.0);
        std::complex <double> sum2(0.0,0.0);
        std::complex <double> p1=un/(un-pow(deux,un-c1));

        for(int k=1;k<=N;k++){
                 std::complex <double> p2=un/pow(k,c1);
                 if(k%2==0)sum1+=p2;
                 if(k%2==1)sum1-=p2;
        }
        std::vector<double   > V1(N);
        std::vector<double   > V2(N);
        double coef=1.0;
        double up=N;
        double dw=1.0;
        double su=0.0;
        for(int k=0;k<N;k++){
                coef*=up;up-=1.0;
                coef/=dw;dw+=1.0;
                V1[k]=coef;
                su+=coef;
        }
        for(int k=0;k<N;k++){
                V2[k]=su;
                su-=V1[k];
        }
        for(int k=N+1;k<=2*N;k++){
                 std::complex <double> p2=un/pow(k,c1);
                 double ek=V2[k-N-1];
                 std::complex <double> c3(ek,0.0);
                 std::complex <double> c4=p2*c3;
                 if(k%2==0)sum2+=c4;
                 if(k%2==1)sum2-=c4;
        }

        std::complex <double> rez=(sum1+sum2/pow(deux,N))*p1;
        return(rez);
}

void test_one_zero(double t)
{
	double RS=Z(t,4);
	std::complex <double> c1=test_zerod(t,10);
	std::complex <double> c2=test_zerod(t,100);
	std::complex <double> c3=test_zerod(t,1000);
	std::cout << std::setprecision(15);
        std::cout << "RS= "<<" "<<RS<<" TEST10= "<< c1 << " TEST100=" << c2 << " TEST1000=" << c3 << std::endl;
	
}

void tests_zeros()
{
	test_one_zero(14.1347251417346937904572519835625);
        test_one_zero(101.3178510057313912287854479402924);
        test_one_zero(1001.3494826377827371221033096531063);
        test_one_zero(10000.0653454145353147502287213889928);

}

/*
	An option to better the performance of Z(t) for large values of t is to simplify the equations
	to validate we present a function that tests the known zeros :  look at https://www.lmfdb.org/zeros/zeta/?limit=10&N=10
	We should obtain 0.0
        no need to test many zeros. In case of a bug the column 2 will show large values instead of values close to 0 like with the original code
	Observe that when t increases the accuracy increases until the limits of the IEEE 754 norm block us, we should work with extended precision
	But here a few digits of precision are enough to count the zeros, only on rare cases the _float128 should be used
	But this limitation only appears very far and with the constraint of resources it won't be possible to reach this region. 
	----------------------------------------------------------------------------------------------------------------------
	value in double			should be 0.0		 value in strings: LMFDB all the digits are corrects
        14.13472514173469463117        -0.00000248590756340983   14.1347251417346937904572519835625
        21.02203963877155601381        -0.00000294582959536882   21.0220396387715549926284795938969
        25.01085758014568938279        -0.00000174024500421144   25.0108575801456887632137909925628
       178.37740777609997167019         0.00000000389177887139   178.3774077760999772858309354141843
       179.91648402025700193008         0.00000000315651035865   179.9164840202569961393400366120511
       182.20707848436646258961         0.00000000214091858131   182.207078484366461915407037226988
 371870901.89642333984375000000         0.00000060389888876036   371870901.8964233245801283081720385309201
 371870902.28132432699203491211        -0.00000083698274928878   371870902.2813243157291041227177012243450
 371870902.52132433652877807617        -0.00000046459056067712   371870902.5213243412580878836297930128983
*/

char line[1024];
void test_fileof_zeros(const char *fname)
{
	FILE *fi=fopen(fname,"r");
	assert(fi!=NULL);
	for(;;){
		double t,RS;
		fgets(line,1000,fi);
		if(feof(fi))break;
		sscanf(line,"%lf",&t);
		RS=Z(t,4);
		printf(" %30.20lf %30.20lf   %s",t,RS,line);

	}
	fclose(fi);
}

int main(int argc,char **argv)
{
	double LOWER,UPPER,SAMP;
	const double pi = 3.1415926535897932385;
	//tests_zeros();
	//test_fileof_zeros("ZEROS");
	try {
		LOWER=std::atof(argv[1]);
		UPPER=std::atof(argv[2]);
		SAMP =std::atof(argv[3]);
	}
	catch (...) {
                std::cout << argv[0] << " START END SAMPLING" << std::endl;
                return -1;
        }
	double estimate_zeros=theta(UPPER)/pi;
	printf("I estimate I will find %1.3lf zeros\n",estimate_zeros);

	double STEP = 1.0/SAMP;
	ui64   NUMSAMPLES=floor((UPPER-LOWER)*SAMP+1.0);
	double prev=0.0;
	double count=0.0;
	double t1=dml_micros();
	for (double t=LOWER;t<=UPPER;t+=STEP){
		double zout=Z(t,4);
		if(t>LOWER){
			if(   ((zout<0.0)and(prev>0.0))
			    or((zout>0.0)and(prev<0.0))){
				//printf("%20.6lf  %20.12lf %20.12lf\n",t,prev,zout);
				count++;
			}
		}
		prev=zout;
	}
	double t2=dml_micros();
	printf("I found %1.0lf Zeros in %.3lf seconds\n",count,(t2-t1)/1000000.0);
	return(0);
}



