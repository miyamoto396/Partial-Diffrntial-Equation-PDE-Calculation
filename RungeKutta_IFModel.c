/***************************************
RungeKutta_IFModel
20151107
ver0.0
***************************************/
/*--------------------------------------
脳の計算論p6より
Integrate-and-fire(積分発火型)モデル
の計算を前fortranでつくったルンゲクッタ法で計算してみる
----------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*初期値*/
#define tini 0.0
#define Vini 0.0		/*膜の初期電位*/
#define tfin 10.0
#define tNbin 1000		/*刻み数*/
#define tau 1.0			/*膜の時定数*/
#define Etau 0.0		/*たぶん膜外の電位*/

double function();

/*定数を定義*/
double pi=M_PI;		/*math.hの中のM_PI*/
double e=M_E;		/*自然対数*/

int main()
{
	FILE *fp;						/*ファイルポインタの宣言*/
	double t,V,I,dt,VMax;
	double k1,k2,k3,k4;
	int i;

	VMax=0.5;						/*閾値*/
	dt=(tfin-tini)/tNbin;			/*刻み幅*/
	t=tini;							/*時間の初期値を代入*/
	V=Vini;							/*膜電圧の初期値を代入*/

	/*ファイルオープン*/
	if ((fp = fopen("IFModel.d", "w")) == NULL) {
		printf("file open error!!\n");
		exit(EXIT_FAILURE);
	}

	/*ルンゲクッタ法による計算*/
	for(i=1;i<=tNbin;i++){

		fprintf(fp, "%f %f\n",t,V);

		k1=dt*function(t,V);
		k2=dt*function(t+0.5*dt,V+0.5*k1);
		k3=dt*function(t+0.5*dt,V+0.5*k2);
		k4=dt*function(t+dt,V+k3);
		V=V+(k1+2.0*k2+2.0*k3+k4)/6.0;
		t=tini+dt*i;
		
		/*閾値を超えたら初期値をいれる*/
		if(V>VMax){
			V=Vini;
		}
	}

	fclose(fp);				/*ファイルクローズ*/
	return 0;
}

/*積分発火型(IF)モデル*/
double function(double t,double V)
{
	double dVdt,I;

	I=1.0;		/*入力電流*/

	dVdt=(-V+Etau+I)/tau;

	return dVdt;
}