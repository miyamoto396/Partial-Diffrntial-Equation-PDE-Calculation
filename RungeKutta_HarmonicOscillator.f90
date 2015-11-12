!-------------------
!
!20150820
!
!runge_kutta_method
!Harmonic_Oscillator
!
!funcdy:dy/dx=F1(x,y,z)速度にあたる
!funcdz:dz/dx=F2(x,y,z)解きたい2次の微分方程式
!
!-------------------

program main
implicit none

!set parameter
double precision,parameter::tini=0.0d0		!t初期値(時間)
double precision,parameter::xini=0.0d0		!x初期値(位置)
double precision,parameter::vini=1.0d0		!v初期値(初速度)


double precision,parameter::tfin=10.0d0		!xをどこまで計算するか
!double precision,parameter::zfin=50.0d0
integer,parameter::tNbin=10000				!刻み数

integer::i
double precision::t,x,v
double precision::dt
double precision::k1,k2,k3,k4,q1,q2,q3,q4
open(10,file='RungeKutta_Harmo.d')

dt=(tfin-tini)/dble(tNbin)		!刻む幅を計算
t=tini 							!初期値を代入
x=xini							!初期値を代入
v=vini

!ルンゲクッタ法の計算をする
do i=1,tNbin
	k1=dt*funcdy(t,x,v)
	q1=dt*funcdz(t,x,v)

	k2=dt*funcdy(t+0.50d0*dt,x+0.50d0*k1,v+0.50d0*q1)
	q2=dt*funcdz(t+0.50d0*dt,x+0.50d0*k1,v+0.50d0*q1)
	
	k3=dt*funcdy(t+0.50d0*dt,x+0.50d0*k2,v+0.50d0*q2)
	q3=dt*funcdz(t+0.50d0*dt,x+0.50d0*k2,v+0.50d0*q2)
	
	k4=dt*funcdy(t+dt,x+k3,v+q3)
	q4=dt*funcdz(t+dt,x+k3,v+q3)

	x=x+(k1+2d0*k2+2d0*k3+k4)/6d0
	v=v+(q1+2d0*q2+2d0*q3+q4)/6d0
    t=tini+dt*dble(i)
    write(10,'(2e26.16)') t,x,v
end do

close(10)

!-----------------------------------------

stop
contains


!dy/dx
function funcdy(T,X,V) result(dydx)		
	double precision,intent(in)::T,X,V
	double precision dydx

	dydx=V  		

end function funcdy

!dz/dx
function funcdz(T,X,V) result(dzdx)
	double precision,intent(in)::T,X,V
	double precision dzdx

	dzdx=-X  		

end function funcdz

end program main