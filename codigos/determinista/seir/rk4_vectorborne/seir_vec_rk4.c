// COMPILAÇÃO: gcc seir_vec_rk4.c -Wall -lm -o seir_vec_rk4
// EXECUÇÃO: ./seir_vec_rk4 s0 e0 i0 r0 sv0 ev0 iv0 dt tf np (suscetíveis, infectados e vacinados iniciais, intervalo de tempo, tempo total e número de pontos)
// valores recomendados: s0= 0.999, e0= 0.0, i0=0.001, r0= 0.0, sv0= 0.99, ev0= 0.0, iv0= 0.01, dt= 0.1, tf= 200.0, np= 2000
// 1 = Suscetível, 2 = Exposto, 3 = Infectado, 4 = Recuperado, 5 = Suscetível Vetor, 6 = Exposto Vetor, 7 = Infectado Vetor

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

double f(double S, double Iv){ // primeira eq dif S
	double beta= 0.35;
	double B = 0.02;
	double mu_h = 0.02;
	return (-beta*S*Iv + B - mu_h*S);
}

double g(double S, double Iv, double E){ // segunda eq dif E
	double beta= 0.35;
	double kappa= 0.1;
	double mu_h = 0.02;
	return (beta*S*Iv - kappa*E - mu_h*E);
}

double h(double E, double I){ // terceira eq dif I
	double gamma= 0.14;
	double kappa= 0.1;
	double mu_h = 0.02;
	return (kappa*E - gamma*I - mu_h*I);
}

double fr(double I, double R){ // terceira eq dif R
	double gamma= 0.14;
	double mu_h = 0.02;
	return (gamma*I - mu_h*R);
}

double fsv(double Sv, double I){ // primeira eq dif Sv
	double A= 0.02;
	double mu= 0.02;
	double beta_v= 0.15;
	return (A - mu*Sv - beta_v*Sv*I);
}

double fse(double Sv, double I, double Ev){ // primeira eq dif Sv
	double mu= 0.02;
	double eta= 0.2;
	double beta_v= 0.15;
	return (beta_v*Sv*I - mu*Ev - eta*Ev);
}

double fsi(double Iv, double Ev){ // primeira eq dif Sv
	double mu= 0.02;
	double eta= 0.2;
	return (eta*Ev - mu*Iv);
}

int main(int argc, char **argv){
	if (argc != 11){
		printf("Erro, da uma lida na execução do programa, número errado de entradas.\n");
		exit(1);
	};

	int i, j, k, n, np;
	double S, E, I, R, Sv, Ev, Iv, t, tf, dt, k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4, nn1, nn2, nn3, nn4, kv1, kv2, kv3, kv4, lv1, lv2, lv3, lv4, mv1, mv2, mv3, mv4;

	S= atof(argv[1]);
	E= atof(argv[2]);
	I= atof(argv[3]);
	R= atof(argv[4]);
	Sv= atof(argv[5]);
	Ev= atof(argv[6]);
	Iv= atof(argv[7]);
	dt= atof(argv[8]);
	tf= atof(argv[9]);
	np= atoi(argv[10]);
	n=tf/(np*dt);
	FILE *arq;
	arq= fopen("seir_vec_rk4.dat", "w");
	t= 0.0;

	for(i= 0; i< np; i++){
		for(j= 0; j< n; j++){ 
			k1= f(S, Iv); // inclinação da reta tangente
			l1= g(S, Iv, E); // inclinação da reta tangente
			m1= h(E, I);
			nn1= fr(I, R);
			kv1= fsv(Sv, I);
			lv1= fse(Sv, I, Ev);
			mv1= fsi(Iv, Ev);
			k2= f(S+0.5*dt*k1, Iv+0.5*dt*mv1);
			l2= g(S+0.5*dt*k1, Iv+0.5*dt*mv1, E+0.5*dt*l1);
			m2= h(E+0.5*dt*l1, I+0.5*dt*m1);
			nn2= fr(I+0.5*dt*m1, R+0.5*dt*nn1);
			kv2= fsv(Sv+0.5*dt*kv1, I+0.5*dt*m1);
			lv2= fse(Sv+0.5*dt*kv1, I+0.5*dt*m1, Ev+0.5*dt*lv1);
			mv2= fsi(Iv+0.5*dt*mv1, Ev+0.5*dt*lv1);
			k3= f(S+0.5*dt*k2, Iv+0.5*dt*mv2);
			l3= g(S+0.5*dt*k2, Iv+0.5*dt*mv2, E+0.5*dt*l2);
			m3= h(E+0.5*dt*l2, I+0.5*dt*m2);
			nn3= fr(I+0.5*dt*m2, R+0.5*dt*nn2);
			kv3= fsv(Sv+0.5*dt*kv2, I+0.5*dt*m2);
			lv3= fse(Sv+0.5*dt*kv2, I+0.5*dt*m2, Ev+0.5*dt*lv2);
			mv3= fsi(Iv+0.5*dt*mv2, Ev+0.5*dt*lv2);
			k4= f(S+dt*k3, Iv+dt*mv3);
			l4= g(S+dt*k3, Iv+dt*mv3, E+dt*l3);
			m4= h(E+dt*l3, I+dt*m3);
			nn4= fr(I+dt*m3, R+dt*nn3);
			kv4= fsv(Sv+dt*kv3, I+dt*m3);
			lv4= fse(Sv+dt*kv3, I+dt*m3, Ev+dt*lv3);
			mv4= fsi(Iv+dt*mv3, Ev+dt*lv3);
			S+= dt*(k1+2.0*k2+2.0*k3+k4)/6.0; //passo temporal
			E+= dt*(l1+2.0*l2+2.0*l3+l4)/6.0; //passo temporal
			I+= dt*(m1+2.0*m2+2.0*m3+m4)/6.0;
			R+= dt*(nn1+2.0*nn2+2.0*nn3+nn4)/6.0;
			Sv+= dt*(kv1+2.0*kv2+2.0*kv3+kv4)/6.0;
			Ev+= dt*(lv1+2.0*lv2+2.0*lv3+lv4)/6.0;
			Iv+= dt*(mv1+2.0*mv2+2.0*mv3+mv4)/6.0;
			t+= dt;
		}
		fprintf(arq, "%e %e %e %e %e %e %e %e\n", t, S, E, I, R, Sv, Ev, Iv);
	}

	fclose(arq);
	return 0;
}
