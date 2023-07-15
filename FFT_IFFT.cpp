#include <stdio.h>
#include <math.h>

typedef struct complex
{
    double real;
    double imag;
}complex;

void read_txtFile(int N, complex x[], char fn[]){
    FILE *fp;
    int i;
    if((fp=fopen(fn,"r"))==NULL){
        printf("failed\n");
    }else{
        for(i=0;i<N;i++){
            fscanf(fp,"%lf",&x[i].real);
        }
    }
    fclose(fp);
}

void write_txtFile(int N, complex F[], char fn[]){
    FILE *fp;
    int i;
    if((fp=fopen(fn,"w"))==NULL){
        printf("write process is failed.\n");
    }else{
        for(i = 0; i<N ;i++){
            fprintf(fp,"%lf %lf\n",F[i].real,F[i].imag);
        }
    }
}


complex add_complex(complex a, complex b){
    complex re_complex;
    re_complex.real = a.real + b.real;
    re_complex.imag = a.imag + b.imag;
    return re_complex;
}

complex sub_complex(complex a, complex b){
    complex re_complex;
    re_complex.real = a.real - b.real;
    re_complex.imag = a.imag - b.imag;
    return re_complex;
}

complex mul_complex(complex a, complex b){
    complex re_complex;
    re_complex.real = a.real*b.real - a.imag*b.imag;
    re_complex.imag = a.real*b.imag + a.imag*b.real;
    return re_complex;
}
// Twiddle factor for DFT
complex W(int N, int x){
    complex re_complex;
    re_complex.real = cos(2*M_PI*x/N);
    re_complex.imag = -sin(2*M_PI*x/N);
    return re_complex;
}
// Twiddle factor for IDFT
complex W_i(int N, int x){
    complex re_complex;
    re_complex.real = cos(2*M_PI*x/N);
    re_complex.imag = sin(2*M_PI*x/N);
    return re_complex;
}

int Bitreverse(int N,int i){
    int r_i = 0x00;

    int cnt = N;
    int flag = 0x01|cnt;
    int mask = 0x1;
    int a;
    while (flag!=1)
    {
        cnt = cnt >> 1;
        flag = 0x01|cnt;

        a = i&cnt;
        if(a != 0x00){
            r_i = (r_i | mask); 
        }

        mask = mask << 1;
    }
    return r_i;
}
// dft_idft=>1:DFT, 0:IDFT
void FFT_Butterfly_Recursion(int N, complex x[], complex X[],int dft_idft){

    int i;
    complex x_next_upper[N/2];
    complex x_next_below[N/2];

    // complex X_next_upper[N/2];
    // complex X_next_below[N/2];
    complex *X_next_upper;
    complex *X_next_below;
    X_next_upper = X;
    X_next_below = &X[N/2];

    for ( i = 0; i < N/2; i++)
    {
        x_next_upper[i] = add_complex(x[i],x[i+N/2]);
        switch (dft_idft)
        {
        case 1:
            x_next_below[i] = mul_complex(W(N,i),sub_complex(x[i],x[i+N/2]));    
            break;
        case 0:
            x_next_below[i] = mul_complex(W_i(N,i),sub_complex(x[i],x[i+N/2])); 
            break;   
        default:
            x_next_below[i] = mul_complex(W(N,i),sub_complex(x[i],x[i+N/2]));    
            break;
        }
    }

    if(N==2){
        X[0] = x_next_upper[0];
        X[1] = x_next_below[0];
        return;
    }

    FFT_Butterfly_Recursion(N/2,x_next_upper,X_next_upper,dft_idft);
    FFT_Butterfly_Recursion(N/2,x_next_below,X_next_below,dft_idft);

    return;

}


void FFT(int N, complex x[], complex X[]){
    int i;
    complex X_tmp[N];

    FFT_Butterfly_Recursion(N,x,X_tmp,1);

    for ( i = 0; i < N; i++)
    {
        X[i] = X_tmp[Bitreverse(N,i)];
    }
}

void IFFT(int N, complex x[], complex X[]){
    int i;
    complex X_tmp[N];

    FFT_Butterfly_Recursion(N,x,X_tmp,0);

    for ( i = 0; i < N; i++)
    {
       X[i].real = X_tmp[Bitreverse(N,i)].real/N;
       X[i].imag = X_tmp[Bitreverse(N,i)].imag/N;
    }
}

int main(int argc, char const *argv[])
{
    int N = 256;
    complex input_data[N];
    complex result[N];
    complex result_r[N];


    char read_fileName[]   = "readText_sample.txt";
    char write_fileName[]  = "writeText_sample_FFT.txt";
    char write_fileName2[] = "writeText_sample_IFFT.txt";


    // input data generation
    // int i;
    // for(i = 0; i<N;i++){
    //     data[i].real = 30 * sin(5*(2*M_PI/N)*i)+5*cos(30*(2*M_PI/N)*i);
    //     data[i].imag = 0.0;
    // }
    read_txtFile(N,input_data,read_fileName);

    FFT(N,input_data,result);     

    IFFT(N,result,result_r);

    write_txtFile(N,result,write_fileName);
    write_txtFile(N,result_r,write_fileName2);

    return 0;
}
