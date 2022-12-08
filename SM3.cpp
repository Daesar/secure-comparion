#include "SM3.h"

using namespace std;


char message_buffer[64]={0};
uint32_t  T[64]={0};


// void out_hex()
// {
// 	for (int i = 0; i < 8; i++)
// 	{
// 		cout<<hex<<hashV[i]<<"  ";
// 	}
// 	cout<<endl;
	
// }



//向左循环移k位
uint32_t rotate_left(uint32_t a,uint32_t k)
{
	k = k % 32;
	return ((a << k) & 0xFFFFFFFF) | ((a & 0xFFFFFFFF) >> (32 - k));
}


int init_T(){
    int i=0;
    for(i=0;i<16;i++){
        T[i]=0x79cc4519;
    }
    for(i=16;i<64;i++){
        T[i]=0x7a879d8a;
    }
    return 1;
}



uint32_t FF(uint32_t X,uint32_t Y,uint32_t Z,int j){
    if(0<=j&&j<16){
        return X^Y^Z;
    }
    else if(j>=16&&j<64){
        return (X&Y)|(X&Z)|(Y&Z);
    }
	else return -1;
}
uint32_t GG(uint32_t X,uint32_t Y,uint32_t Z,int j){
    if(j>=0&&j<16){
        return X^Y^Z;
    }
    else if(j>=16&&j<64){
        return (X&Y)|((~X)&Z);
    }
	else return -1;
}



uint32_t P0(uint32_t X){
    return (X ^ (rotate_left(X, 9)) ^ (rotate_left(X, 17)));
}
uint32_t P1(uint32_t X){
    return ( X ^ (rotate_left(X, 15)) ^ (rotate_left(X, 23)));
}


int CF(char *arr,uint32_t *hashV){
    uint32_t W[68];
    uint32_t W1[64];
    uint32_t A,B,C,D,E,F,G,H;
    uint32_t SS1,SS2,TT1,TT2;


    //消息拓展
    for(int j=0;j<16;j++){
        W[j]=(arr[j * 4 + 0] << 24) | (arr[j * 4 + 1] << 16) | (arr[j * 4 + 2] << 8) | (arr[j * 4 + 3]);
    }
    for (int j = 16; j < 68; j++)
	{
		W[j] = P1(W[j - 16] ^ W[j - 9] ^ (rotate_left(W[j - 3], 15))) ^ (rotate_left(W[j - 13], 7)) ^ W[j - 6];
	}
	for (int j = 0; j < 64; j++)
	{
		W1[j] = W[j] ^ W[j + 4];
	}

    A=hashV[0];
    B=hashV[1];
    C=hashV[2];
	D=hashV[3];
	E=hashV[4];
	F=hashV[5];
	G=hashV[6];
	H=hashV[7];

    for (int j = 0; j < 64; j++)
	{
		SS1 = rotate_left(((rotate_left(A, 12)) + E + (rotate_left(T[j], j))) & 0xFFFFFFFF, 7);
		SS2 = SS1 ^ (rotate_left(A, 12));
		TT1 = (FF(A, B, C, j) + D + SS2 + W1[j]) & 0xFFFFFFFF;
		TT2 = (GG(E, F, G, j) + H + SS1 + W[j]) & 0xFFFFFFFF;
		D = C;
		C = rotate_left(B, 9);
		B = A;
		A = TT1;
		H = G;
		G = rotate_left(F, 19);
		F = E;
		E = P0(TT2);

	}
    hashV[0] = (A ^ hashV[0]);
	hashV[1] = (B ^ hashV[1]);
	hashV[2] = (C ^ hashV[2]);
	hashV[3] = (D ^ hashV[3]);
	hashV[4] = (E ^ hashV[4]);
	hashV[5] = (F ^ hashV[5]);
	hashV[6] = (G ^ hashV[6]);
	hashV[7] = (H ^ hashV[7]);
	return 1;
}

void SM3_Init(uint32_t *hashV)
{
	
	init_T();
	hashV[0] = 0x7380166f;
	hashV[1] = 0x4914b2b9;
	hashV[2] = 0x172442d7;
	hashV[3] = 0xda8a0600;
	hashV[4] = 0xa96f30bc;
	hashV[5] = 0x163138aa;
	hashV[6] = 0xe38dee4d;
	hashV[7] = 0xb0fb0e4e;
}



void Fill(char *msg,int msglen,uint32_t *hashV){
	int left = 0;
	int i=0;
	uint64_t total = msglen * 8;
		
	for(i = 0; i < msglen/64; i++){
		memcpy(message_buffer, msg + i * 64, 64);
		CF(message_buffer,hashV);
	}


	left = msglen%64;
	memcpy(message_buffer, msg + i * 64, left);
	//填充0
	memset(&message_buffer[left], 0, 64 - left);	
	message_buffer[left] = 0b10000000;
	if(left <= 55){
		for (i = 0; i < 8; i++)
			message_buffer[56 + i] = (total >> ((7 - i) * 8)) & 0xFF;
		CF(message_buffer,hashV);
	}else{
		CF(message_buffer,hashV);
		memset(message_buffer, 0, 64);
		for (i = 0; i < 8; i++)
			message_buffer[56 + i] = (total >> ((8 - i) * 8)) & 0xFF;
		CF(message_buffer,hashV);
	}
}

uint32_t *SM3( char *msg,int msglen)
{
	uint32_t *hashV=(uint32_t *)malloc(8*sizeof(uint32_t));
	SM3_Init(hashV);
	Fill(msg,msglen,hashV);
	//out_hex();
	return hashV;
}
// int main(){
// 	//test1: msg1="abc"
// 	cout<<"meg1=abc"<<endl;
// 	unsigned char msg1[]="abc";
// 	SM3(msg1,3);

// 	//test2: meg2=61626364 61626364 61626364 61626364 61626364 61626364 61626364 61626364
// 	//            61626364 61626364 61626364 61626364 61626364 61626364 61626364 61626364
// 	cout<<"meg2=61626364 61626364 61626364 61626364 61626364 61626364 61626364 6162636461626364 61626364 61626364 61626364 61626364 61626364 61626364 61626364"<<endl;
// 	unsigned char msg2[]={0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 
// 	                      0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64,
// 						  0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64,
// 						  0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64, 0x61,0x62,0x63,0x64};
// 	SM3(msg2,64);
// }