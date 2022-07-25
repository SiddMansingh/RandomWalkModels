#include <stdio.h>
 #include <math.h>
#include<iostream>
#include <sstream>
#include<fstream>
#include <cstdlib>
#include <iomanip>
#define PI 3.141592653589793238462643383279502884L
using namespace std; 
int l1 = 0;
int l2 = 0;
long long lA=0;
long long lB=0;
long long cons=0;
int n1,n2,n,k,iterations;
long long tw,size,size0,size1;
int ch;
double j2 ;
long long tempo;
long long count;
int vals=10;
int omegaiters = 2001;
int lim=4;
double epsilon = 0.005;
double E0;
long long int* lookupA0 = NULL;
long long int* lookupB0 = NULL;
long long int* lookupA1 = NULL;
long long int* lookupB1 = NULL;
long long int* states0 = NULL;
long long int* states1 = NULL;
int** cmatrix = NULL;
double* cstrength = NULL;
double* psi0 = NULL;
double** psi1 = NULL;
double** phi0 = NULL;
double** phi1 = NULL;
double* Lancz_a = NULL;
double* Lancz_b = NULL;
double* scatter = NULL;
double** kappa = NULL;
double normphi0;
ofstream scat_file;


long long int binomialCoeff(int n, int k){                                              //Computing Binomial Coefficients
  if (k==0 || k==n)
    return 1;
  return  binomialCoeff(n-1, k-1) + binomialCoeff(n-1, k);
}
unsigned int countSetBits(unsigned int n){                                              //Count number of Set Bits (Useful for checking which space the state is in
  unsigned int count = 0;
  while (n)
  {
    count += n & 1;
    n >>= 1;
  }
  return count;
}
void lookup_create_A(int a[],int l1,int l2,long long lookupA[],long long lookupB[] ){
	int i=0;
	lA=0;
    lB=0;
    while(a[i]<n/2 && i<k){
    	lA=lA+pow(2,n/2-1-a[i]);
    	i++;
	}
	while(i<k){
		lB=lB+pow(2,(n-a[i]-1));
		i++;
	}
	
	lookupA[lA]=l1;
	lookupB[lB]=l2;
}    
void lookup_create_B(int a[],int l2,long long lookupB[]){
	int i=k-1;
	lB=0;
	while(a[i]>=n/2){
		lB=lB+pow(2,(n-a[i]-1));
		i--;
	}
	lookupB[lB]=l2;
}
void print_fn(int a[], int l1, int l2, long long lookupA[], long long lookupB[]){
	string Dec2Bin, Dec2BinA, Dec2BinB;
	lA = 0;
	lB = 0;
	long long lnet = 0;
	for(int i=0;i<n;i++)
		Dec2Bin += '0';
	for(int i=0;i<n;i++){
		if(i<k)
		Dec2Bin[a[i]] = '1';
		lnet += pow(2,a[i]);
		if(a[i]<n/2 && i<k){
			lA=lA+pow(2,n/2-1-a[i]);
			//Dec2BinA += to_string(a[i]);
		}
			
		else
			if(i<k){
				//Dec2BinB += to_string(a[i]);
				lB=lB+pow(2,(n-a[i]-1));
			}
			
	}
	cout<<lnet<<"\t\t"<<Dec2Bin<<"\t\t"<<Dec2Bin.substr(0,n/2)<<"\t\t"<<lA<<"\t\t"<<lookupA[lA]<<"\t\t"<<Dec2Bin.substr(n/2)<<"\t\t"<<lB<<"\t\t"<<lookupB[lB]<<
	"\t\t"<<lookupA[lA]+lookupB[lB]<<endl;
}
int next_config(int a[],long long lookupA[],long long lookupB[]) {                      //Generates next Configuration
    int i = k - 1;
    int l;
    long long f;
    a[i]++;
    while ((i >= 0) && (a[i] >= n - k + 1 + i)) {
        --i;
        a[i]++;
    }
 	
    if (a[0] > ch ){										//when the first most significant bit cannot be accomodated in the n-bit string, time to stop
    	l1+=l2+1;
    	return 0;
	}	
    l=i;													//signifies the last updated location in the bit string
    for (i = i + 1; i< k; ++i)
        a[i] = a[i - 1] + 1;
    
    if(a[l]<=n/2){											//if the last updated location was a part of the MSB, all possible configurations have been explored 
    														//in the LSB
    	l1+=l2+1;
    	l2=0;
    	lookup_create_A(a,l1,l2,lookupA,lookupB);
	}
	else{
    	l2++;
    	lookup_create_B(a,l2,lookupB);
	}	
	print_fn(a,l1,l2,lookupA,lookupB);
    return 1;
}
void basis_create(long long * &lookupA,long long * &lookupB,long long * &states){       //Basis Generation
	l1=0;
	l2=0;
	size=pow(2ULL,n/2);									//lookup tables are split into two halves, hence the n/2
	tw=pow(2ULL,n);										//total possible configurations with n spins, irrespective of the symmetry
	lookupA= new long long [size];
	lookupB= new long long [size];
	int a[n]; 											//Essential to creating all possible combinations, carries the location of k up spins
	for (int i=0;i<size;i++)							//Assigning a default value, has no utility
		lookupA[i]=i;
	
	size=binomialCoeff(n,k); 							//For the Sz symmetry, this variable represents the number of states with k up-spins
    if (k!=n/2){
	ch=n-k;												//Comes in handy when trying to restrict generation of combinations
	} 
	else {
		int ch=0;
		size=size/2;									//k up followed by n-k down is a different configuration from n-k up followed by k down
														//unless k = n/2 
	} 
	states= new long long [size];
	for (int i = 0; i < k; ++i)							//creates the first state, k consecutive up spins
        a[i] = i;
	lookup_create_A(a,l1,l2,lookupA,lookupB);			//call lookup create after each state creation
	
	tempo=lA*pow(2,n/2)+lB;								//variable to refer to the configuration's decimal representation in 2^n space
	if(tempo>tw/2)
		tempo=tw-1-tempo;								//With Sz symmetry, we are only storing the configuration with the lower magnetisation
	states[0]=tempo;
	count = 1;	
	cout<<"Number"<<"\t\t"<<"Configuration"<<"\t\t"<<"MSB"<<"\t\t"<<"MSB Bin"<<"\t\t"<<"lA"<<"\t\t"<<"LSB"<<"\t\t"<<"LSB Bin"<<"\t\t"<<"lB"<<"\t\t"<<"lA+lB"<<endl;
	print_fn(a,l1,l2,lookupA,lookupB);
    while (next_config(a,lookupA,lookupB)){
    	tempo=lA*pow(2,n/2)+lB;
		if(tempo>tw/2)
			tempo=tw-1-tempo;
		states[count]=tempo;
		count++;	
    }  
}
int main(int argc, char *argv[]) {	
	iterations=10;
	n1=atoi(argv[1]);
	n2=atoi(argv[2]);							//lattice dimensions n1 X n2
	n=n1*n2;
	k=atoi(argv[3]);							//Sz=k
	basis_create(lookupA0,lookupB0,states0);	//Creates basis for Sz=k
	/*size0=size;
	psi0 = new double[size0];
	k=k+1;
	basis_create(lookupA1,lookupB1,states1); 	//Creates basis for Sz=k+1
	size1=size;
	psi_initialiser();
	phi_initialiser();
	create_kappa3();
	Iterate_j2();*/
	return 0;
}
