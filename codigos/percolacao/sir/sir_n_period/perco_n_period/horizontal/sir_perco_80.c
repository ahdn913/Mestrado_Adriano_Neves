// gcc sir_perco_50_teste.c -O3 -Wall -lgsl -lgslcblas -lm -o sir_perco_50_teste
// argv[1] = seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 80 // L = Ni = Nj
#define Ni 80   // tamanho da rede em i
#define Nj 80   // tamanho da rede em j
#define repeticoes 10000

int phi[L*L];

void op(int t){
	int i, j;
	FILE *arq;
	char nome[100];

	sprintf(nome, "dat/sir-%d.dat", t);
	arq= fopen(nome, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq, "\n");
	}
    fclose(arq);
}

void ic(int *phi){ 
	int i, j, cont_inf;

	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
            if(i == Ni/2 && j == Nj/2){
                phi[i*Nj+j]= 2;
				cont_inf= 1;
            }
            else{
                phi[i*Nj+j]= 1;
            }
		}
	}
}

void recur(int j, int i, int v) {
	if(i >= 0 && j >= 0 && i < L && j < L && phi[j*L+i] == 3){ 
		phi[j*L+i] = v;
		recur((j+1), i, v);
		recur((j-1), i, v);
		recur(j, (i+1), v);
		recur(j, (i-1), v);
	}
}

int count_clusters(void){
	int i, j, cls;
	cls= 6400; 
	for(j= 0; j< L; j++){
		for(i= 0; i< L; i++){
			if(phi[j*L+i] != 3) 
				continue ;
			recur(j, i, ++cls);
		}
	}
	return cls;
}

int main(int argc, char **argv){
	int i, j, x; 
	int ativo, passivo, vizinho, pr_int, cont_inf;
	int max, indicador, vizinho_esqu, vizinho_dire, vizinho_cima, vizinho_baix;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, lambda;
	FILE *file;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(lambda= 0.1005; lambda< 0.2425; lambda+= 0.002){
		ic(phi);
		cont_inf= 1;
		op(0);
		pr_int= lambda*10000;
		while(cont_inf > 0){
			for(x= 0; x< Ni*Nj; x++){
				i= gsl_rng_uniform(w)*Ni;
				j= gsl_rng_uniform(w)*Nj;
				ativo= i*Nj+j;
				vizinho_esqu= vizinho_dire= vizinho_cima= vizinho_baix= 1;

				if((j-1) < 0){
					vizinho_esqu= 0;
				}
				if((j+1) >= Nj){
					vizinho_dire= 0;
				}
				if((i-1) < 0){
					vizinho_cima= 0;
				}
				if((i+1) >= Ni){
					vizinho_baix= 0;
				}

				indicador= 0;

				while(indicador != 1){
					vizinho= gsl_rng_uniform(w)*4;
					switch(vizinho){
						case 0:
							if(vizinho_esqu == 1){
								passivo= i*Nj+(j-1); // esqu	
								indicador= 1;
							}
							else{
								indicador= 0;
							}
						break;
						case 1:
							if(vizinho_dire == 1){
								passivo= i*Nj+(j+1); // dire	
								indicador= 1;
							}
							else{
								indicador= 0;
							}
						break;
						case 2:
							if(vizinho_cima == 1){
								passivo= (i-1)*Nj+j; // esqu	
								indicador= 1;
							}
							else{
								indicador= 0;
							}
						break;
						default:
							if(vizinho_baix == 1){
								passivo= (i+1)*Nj+j; // esqu	
								indicador= 1;
							}
							else{
								indicador= 0;
							}
						break;
					}
				}
				acao= gsl_rng_uniform(w);
				if(acao > lambda){
					//infecção
					if(phi[passivo] == 1 && phi[ativo] == 2){
						phi[passivo]= phi[ativo];
					}
				}
				else{
					//recuperação
					if(phi[ativo] == 2){
						phi[ativo]= 3;
					}
				}
			}
			cont_inf= 0;
			for(i= 0; i < Ni; i++){
				for(j= 0; j < Nj; j++){
					if(phi[i*Nj+j]==2){
						cont_inf++;
					}
				}
			}
		}
		op(1);
	//				Percolação
		key= 0;
		max= count_clusters();
		for(i= 6401; i<= max; i++){ 
			a_ih= a_fh= 0;
			for(j= 0; j< L; j++){
				if(phi[j] == i){ 
					a_ih++; 
				}
				if(phi[(L-1)*L+j] == i){ 
					a_fh++; 
				}
			}

			if(a_ih > 0 && a_fh > 0){
				key= 1;
			}
		}
		
		sprintf(name, "dat/%d_%d/percolacao_%d.dat", repeticoes, L, pr_int);
		if(!(file= fopen(name, "a"))){
			printf("cannot open file teste\n");
			exit(0);
		}
		if(key == 1){
			fprintf(file, "1\n");
			printf("percolou\n");
		}
		else{
			fprintf(file, "0\n");
			printf("não percolou\n");
		}
		fclose(file);
	}
	gsl_rng_free(w);
	return 0;
}
