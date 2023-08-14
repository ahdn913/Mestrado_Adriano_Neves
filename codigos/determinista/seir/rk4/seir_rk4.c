// COMPILAÇÃO: gcc seir_rk4.c -Wall -lm -o seir_rk4
// EXECUÇÃO: ./seir_rk4 s0 e0 i0 r0 dt tf np (suscetíveis, infectados e vacinados iniciais, intervalo de tempo, tempo total e número de pontos)
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

double f(double S, double I){ // primeira eq dif S
	double beta= 0.85;
	return (-beta*S*I);
}

double g(double S, double I, double E){ // segunda eq dif E
	double beta= 0.85;
	double kappa= 0.2;
	return (beta*S*I - kappa*E);
}

double h(double E, double I){ // terceira eq dif I
	double kappa= 0.2;
	double gamma= 0.15;
	return (kappa*E - gamma*I);
}

double fr(double I){ // terceira eq dif R
	double gamma= 0.15;
	return (gamma*I);
}

int main(int argc, char **argv){
	if (argc != 8){
		printf("Erro, dá uma lida na execução do programa, número errado de entradas.\n");
		exit(1);
	};

	int i, j, k, n, np;
	double S, E, I, R, t, tf, dt, k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4, nn1, nn2, nn3, nn4;

	S= atof(argv[1]);
	E= atof(argv[2]);
	I= atof(argv[3]);
	R= atof(argv[4]);
	dt= atof(argv[5]);
	tf= atoi(argv[6]);
	np= atoi(argv[7]);
	n=tf/(np*dt);
	FILE *arq;
	arq= fopen("seir_rk4.dat", "w");
	t= 0.0;

	for(i= 0; i< np; i++){
		for(j= 0; j< n; j++){ 
			k1= f(S, I); // inclinação da reta tangente
			l1= g(S, I, E); // inclinação da reta tangente
			m1= h(E, I);
			nn1= fr(I);
			k2= f(S+0.5*dt*k1, I+0.5*dt*l1);
			l2= g(S+0.5*dt*k1, I+0.5*dt*l1, E+0.5*dt*m1);
			m2= h(E+0.5*dt*m1, I+0.5*dt*l1);
			nn2= fr(I+0.5*dt*l1);
			k3= f(S+0.5*dt*k2, I+0.5*dt*l2);
			l3= g(S+0.5*dt*k2, I+0.5*dt*l2, E+0.5*dt*m2);
			m3= h(E+0.5*dt*m2, I+0.5*dt*l2);
			nn3= fr(I+0.5*dt*l2);
			k4= f(S+dt*k3, I+dt*l3);
			l4= g(S+dt*k3, I+dt*l3, E+dt*m3);
			m4= h(E+dt*m3, I+dt*l3);
			nn4= fr(I+dt*l3);
			S+= dt*(k1+2.0*k2+2.0*k3+k4)/6.0; //passo temporal
			E+= dt*(l1+2.0*l2+2.0*l3+l4)/6.0; //passo temporal
			I+= dt*(m1+2.0*m2+2.0*m3+m4)/6.0; //passo temporal
			R+= dt*(nn1+2.0*nn2+2.0*nn3+nn4)/6.0;
			t+= dt;
		}
		fprintf(arq, "%e %e %e %e %e\n", t, S, E, I, R);
	}

	fclose(arq);
	return 0;
}
