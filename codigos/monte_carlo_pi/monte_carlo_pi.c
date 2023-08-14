// Modelo para calcular pi
// compilação: gcc monte_carlo_pi.c -O3 -Wall -lgsl -lgslcblas -lm -o monte_carlo_pi
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define N 20000  

struct intble {
	int    c; 
	double x; 
	double y;
};

int main(int argc, char **argv){
	int i, j, n, cont_d;
	double d, pi, cont_D, NN;
	struct intble *p;
	FILE *arq;
	char clus[100];

	p= (struct intble *) calloc(N*N, sizeof(struct intble)); // Calloc aloca uma matriz na memória (número, tamanho)

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); // inicialização do gerador de número aleatório (algorítmo Taus)
	for(i= 0; i< N; i++){					
		p[i].x= gsl_rng_uniform(w);
		p[i].y= gsl_rng_uniform(w);
	}

	cont_d= 0;
	sprintf(clus, "dat/p%d.dat", N); 
	arq= fopen(clus, "w"); 	
	for(i= 0; i< N; i++){
		d= sqrt(p[i].x*p[i].x + p[i].y*p[i].y);
		if(d <= 1.0){
			cont_d++;
			p[i].c= 1;
		}
		else{
			p[i].c= 2;
		}
		fprintf(arq, "%e %e %d %d %d\n", p[i].x, p[i].y, p[i].c, N, cont_d); 
	}
	cont_D= cont_d;
	NN= N;
	pi= cont_D/NN*4;
	printf("%d %d\n", cont_d, N);
	printf("pi é %e\n", pi);

	gsl_rng_free(w);
	free(p);
	return 0;
}
