#include <iostream>     
#include <gmp.h>       
#include <gmpxx.h>        
#include <chrono>
#include <math.h>
#include <string>
#include <iomanip>
#include <ctime>
#include <stdlib.h>

#include "SM3.h"

using namespace std;



typedef struct RhysPrivateKey
{
    mpz_t x;
}RhysPrivateKey;
typedef  struct  RhysPublicKey
{
    mpz_t n;
    int b;
    int d;
    mpz_t g;
    mpz_t h;
    int u;
}RhysPublicKey;

RhysPrivateKey *InitRhysPrivateKey(){

    RhysPrivateKey *rhysPrivateKey;
    rhysPrivateKey = (RhysPrivateKey *)malloc(sizeof(*rhysPrivateKey));
    return rhysPrivateKey;
}
RhysPublicKey *InitRhysPublicKey(){
    RhysPublicKey *rhysPublicKey;
    rhysPublicKey = (RhysPublicKey *)malloc(sizeof(*rhysPublicKey));
    return rhysPublicKey;
}


void keyGen(RhysPublicKey *rhysPublicKey,RhysPrivateKey *rhysPrivateKey,int pslen,int qslen){//qslen,pslen=512
                                                                                             //u=512
                                                                                             //b=2
                                                                                             //d=600
                                                                                             //N=p*q
                                                                                             //p=2*b^d*ps+1;q=2*b^d*qs+1
                                                                                             //security parameter lambda=256
                                                                                             //a>=lambda=256;a=256
                                                                                             
    mpz_t  p,q,ps,qs,g,h;
	mpz_inits(p,q,ps,qs,g,h,NULL);
    
    gmp_randstate_t grt;
    clock_t time = clock(); 
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, time);

	mpz_urandomb(ps, grt, pslen);
	mpz_urandomb(qs, grt, qslen);
	mpz_nextprime(ps, ps);
    mpz_nextprime(qs, qs); 

    mpz_t temp;
    mpz_init(temp);
    mpz_ui_pow_ui(temp,2,600);
    mpz_mul_ui(temp,temp,2);  //temp=2*b^d

    mpz_mul(p,ps,temp);
    mpz_mul(q,qs,temp);
    mpz_add_ui(p,p,1);
    mpz_add_ui(q,q,1);


	while(mpz_probab_prime_p(p, 25)<1){
        mpz_nextprime(ps,ps);
        mpz_mul(p,ps,temp);
		mpz_add_ui(p,p,1);
	}
    while(mpz_probab_prime_p(q, 25) < 1){
        mpz_nextprime(qs,qs);
        mpz_mul(q,qs,temp);
		mpz_add_ui(q,q,1);
	}

    mpz_t n;
    mpz_init(n);
    mpz_mul(n,p,q);    // n=p*q  


    mpz_t x,y,p1,q1,a,div,gp,gq,bd;
    mpz_inits(x,y,p1,q1,a,div,gp,gq,bd,NULL);
   
    mpz_sub_ui(p1,p,1);

    //求gp
	mpz_urandomm(x, grt,p);
    while (!(mpz_cmp_ui(x,1))||!(mpz_cmp_ui(x,0))||!(mpz_cmp(x,p1)))
    {
        mpz_urandomm(x,grt,p);
    }

    mpz_cdiv_q_ui(div,p1,2);
    mpz_powm(y,x,div,p);

    while(!mpz_cmp_ui(y, 1)){
		mpz_add_ui(x, x, 2);
		mpz_powm(y, x, div, p);
	}
    mpz_ui_pow_ui(bd,2,600);
    mpz_cdiv_q(gp,p1,bd);
 
    mpz_powm(gp,x,gp,p); 

    //求在Zq上的gq

    mpz_urandomm(x, grt,q);
    mpz_sub_ui(q1,q,1);
    while (!(mpz_cmp_ui(x,1))||!(mpz_cmp_ui(x,0))||!(mpz_cmp(x,q1)))
    {
        mpz_urandomm(x,grt,q);
    }

    mpz_cdiv_q_ui(div,q1,2);
    mpz_powm(y,x,div,q);

    while(!mpz_cmp_ui(y, 1)){
		mpz_add_ui(x, x, 2);
		mpz_powm(y, x, div, q);
	}

    mpz_cdiv_q(gq,q1,bd);
    mpz_powm(gq,x,gq,q);  

    //用CRT求g
    mpz_t p_Invert,q_Invert;
    mpz_inits(p_Invert,q_Invert,NULL);
    mpz_invert(p_Invert,p,q);
    mpz_invert(q_Invert,q,p);

    mpz_mul(gp,gp,q);
    mpz_mul(gp,gp,q_Invert);
    mpz_mod(gp,gp,n);
    mpz_mul(gq,gq,p);
    mpz_mul(gq,gq,p_Invert);
    mpz_mod(gq,gq,n);

    mpz_add(g,gp,gq);
    mpz_mod(g,g,n);

    //求h

    mpz_t hps,hqs;
    mpz_inits(hps,hqs,NULL);
    
    

    mpz_urandomm(x,grt,p);
    mpz_powm(y, x, temp,p);
	while(!mpz_cmp_ui(y, 1)){
		mpz_add_ui(x, x, 2);
		mpz_powm(y, x, temp, p);
	}
    mpz_set(hps, y);

    mpz_urandomm(x,grt,q);
    mpz_powm(y, x, temp,q);
	while(!mpz_cmp_ui(y, 1)){
		mpz_add_ui(x, x, 2);
		mpz_powm(y, x, temp, q);
	}
    mpz_set(hqs, y);


//CRT

    mpz_mul(hps,hps,q);
    mpz_mul(hps,hps,q_Invert);
    mpz_mod(hps,hps,n);

    mpz_mul(hqs,hqs,p);
    mpz_mul(hqs,hqs,p_Invert);
    mpz_mod(hqs,hqs,n);

    mpz_add(h,hps,hqs);
    mpz_mod(h,h,n);

    mpz_t px,p2,p3;
    mpz_inits(px,p2,p3,NULL);
    mpz_mul(p3,qs,ps);
    mpz_invert(p2,p3,bd);
    mpz_mul(px,p2,p3);

//public key
    mpz_inits(rhysPublicKey->n,rhysPublicKey->g,rhysPublicKey->h,rhysPrivateKey->x,NULL);
    mpz_set(rhysPublicKey->n,n);
    rhysPublicKey->b=2;
    rhysPublicKey->d=600;
    mpz_set(rhysPublicKey->g,g);
    mpz_set(rhysPublicKey->h,h);
    rhysPublicKey->u=512;

//private key
    mpz_set(rhysPrivateKey->x,px);

}
char *RhysEncryption(RhysPublicKey *rhysPublicKey,int m,mpz_t r){
    mpz_t c;
    mpz_init(c);
    mpz_t n,h,g;
    mpz_inits(n,g,h,NULL);
    mpz_set(h,rhysPublicKey->h);
    mpz_set(g,rhysPublicKey->g);
    mpz_set(n,rhysPublicKey->n);
    int b,d,u;
    b=rhysPublicKey->b;
    d=rhysPublicKey->d;
    u=rhysPublicKey->u;
    

    mpz_t hr,gbm,temp;
    mpz_inits(hr,gbm,temp,NULL);
    mpz_ui_pow_ui(temp,b,m);
    mpz_powm(hr,h,r,n);
    mpz_powm(gbm,g,temp,n);
    mpz_mul(c,hr,gbm);
    mpz_mod(c,c,n);
    char *cipher = mpz_get_str(NULL, 10, c);
    return cipher;  
}

char  *RhysDecryption(RhysPrivateKey *rhysPrivateKey,char* c,mpz_t n){
    mpz_t gbm;
    mpz_init(gbm);

    mpz_t cmpz;
    mpz_init(cmpz);
    mpz_set_str(cmpz, c, 10);
    
    mpz_powm(gbm,cmpz,rhysPrivateKey->x,n);
    char *msg = mpz_get_str(NULL, 10, gbm);
    return msg;
}


int Mess(int *m,int j,int length){//mess(m,j)=m[j+1]*b^(j+1)+m[j+2]*b^(j+2)+...+m[length-1]*b[length-1],where b=2
    int mess=0;
    for(int i=j+1;i<length;i++){
        mess=mess+m[i]*pow(2,i);
    }
    return mess;
}
int Length(int a) {
    int count = 0;
    while (a) {
        count++;
        a = a / 2;
    }
    return count;
}

bool comparison(int m1,int m2,RhysPrivateKey *rhysPrivateKey, RhysPublicKey  *rhysPublicKey){ //m1>m2 return 1,else return 0;
    
    int length=Length(m1)>=Length(m2)?Length(m1):Length(m2);
    mpz_t n,h,g;
    mpz_inits(n,g,h,NULL);
    mpz_set(h,rhysPublicKey->h);
    mpz_set(g,rhysPublicKey->g);
    mpz_set(n,rhysPublicKey->n);
    int b,d,u;
    b=rhysPublicKey->b;
    d=rhysPublicKey->d;
    u=rhysPublicKey->u;

    int *m1Set=(int *) malloc(length*sizeof(int));
    int *m2Set=(int *) malloc(length*sizeof(int));
    for(int i=0;i<length;i++){
        m1Set[i]=m1%2;
        m2Set[i]=m2%2;
        m1=m1/2;
        m2=m2/2;
    }

    mpz_t C[length];
    for(int i=0;i<length;i++){
        mpz_init(C[i]);
    }

    int a=256;    // a is apublic integer is chosen by security parameter

    clock_t time=clock();
    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt,time);

    for(int i=0;i<length;i++){
        char *c;
        mpz_t r1;
        mpz_init(r1);
        mpz_urandomb(r1,grt,u);
        if(!mpz_cmp_ui(r1,0)){
        mpz_urandomb(r1,grt,u);
        }
        c=RhysEncryption(rhysPublicKey,a*m1Set[i],r1); 
        mpz_t c_mpzt,temp,mess_si;
        mpz_inits(c_mpzt,temp,mess_si,NULL);
        mpz_set_str(c_mpzt,c,10); 
        int mess=Mess(m1Set,i,length);
        mess=-mess;
        mpz_set_si(mess_si,mess);
        mpz_powm(temp,g,mess_si,n);
        mpz_mul(C[i],temp,c_mpzt);
        mpz_mod(C[i],C[i],n);

    }
    
    uint32_t D2[length][8];    //存放g^v的hash值
    mpz_t u2,v,r2;
    mpz_inits(u2,v,r2,NULL);
    mpz_t D[length];
    for(int i=0;i<length;i++){
        mpz_init(D[i]);
    }

    for(int i=0;i<length;i++){
        mpz_urandomb(u2,grt,a);
        if(!mpz_cmp_ui(u2,0)){
        mpz_urandomb(u2,grt,a);
        }
        mpz_urandomb(v,grt,d);
        if(!mpz_cmp_ui(v,0)){
        mpz_urandomb(v,grt,d);
        }
        mpz_urandomb(r2,grt,u);
        if(!mpz_cmp_ui(r2,0)){
        mpz_urandomb(u2,grt,u);
        }

        //compute D[i]
        mpz_t temp2,temp3,temp4,temp5,bd;
        mpz_inits(temp2,temp3,temp4,temp5,bd,NULL);
        int mess2=Mess(m2Set,i,length);
        mpz_powm_ui(temp2,g,mess2,n);
       
        mpz_mul(temp2,C[i],temp2);
        mpz_mod(temp2,temp2,n); 
        mpz_ui_pow_ui(temp3,b,d-a*(m2Set[i]+1));
        mpz_mul(temp3,u2,temp3);
        mpz_powm(temp4,g,v,n);
        //hash(g^v)
        char *gv=mpz_get_str(NULL,10,temp4);
        uint32_t *hash=new uint32_t[8];
        hash=SM3(gv,strlen(gv));
        for(int j=0;j<8;j++){
            D2[i][j]=hash[j];
        }
        mpz_powm(temp2,temp2,temp3,n);

        mpz_powm(temp5,h,r2,n);
        mpz_mul(temp4,temp4,temp5);
        mpz_mod(temp4,temp4,n);
        mpz_mul(D[i],temp2,temp4);
        mpz_mod(D[i],D[i],n);
       
    }
    

    //B将(D[i],D2[i])以一个随机全排列发送给A
    int A[length];                       
    for(int i=0;i<length;i++){
        A[i]=rand()%length;
        for(int j=0;j<i;j++){
            if(A[i]==A[j]){
                i--;break;
                }
        }  
    }
    for(int i=0;i<length;i++){
        char *Dmes=mpz_get_str(NULL,10,D[A[i]]);
        char *Cmes=RhysDecryption(rhysPrivateKey,Dmes,n);
        uint32_t *hash=new uint32_t[8];
        hash=SM3(Cmes,strlen(Cmes));
    
        for(int j=0;j<8;j++){
            if((hash[j]==D2[A[i]][j])&&j==7){
                return true;
            }else if(hash[j]==D2[A[i]][j]){
                continue;
            }else{
                break;
            }
        }
    }
    return false;

}

int main(){ 
    //KeyGen
    RhysPrivateKey *rhysPrivateKey=InitRhysPrivateKey();
    RhysPublicKey  *rhysPublicKey=InitRhysPublicKey();
    keyGen(rhysPublicKey,rhysPrivateKey,512,512);  //rhys keygen

    srand((unsigned int)time(NULL));

    bool result =comparison(483647252,452151324,rhysPrivateKey, rhysPublicKey);

    if(result){
        cout<<1;
    }else{  
        cout<<0;
    }
    cout<<"over";

    return 1;

}

