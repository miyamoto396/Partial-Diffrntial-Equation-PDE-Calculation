!***********************************
!放物型偏微分方程式(parabolic partial differential equation)を
!差分法(finite difference methods)の陽解法で求める
!20151005
!ver0.0
!***********************************
!(やっぱり陽解法だと刻み幅に制限があってめんどくさい)
!
program PPDE_FDM
implicit none
!--------------------変数定義------------------------
integer,parameter::Nt=1000							!時間刻み数
integer,parameter::Nx=10							!x方向刻み数
integer,parameter::size=1 							!長さ
integer,parameter::time=10 							!時間
double precision,parameter::TempCond=0.07d0			!温度伝導率TemperatureConductivity
double precision,dimension(:,:),allocatable::U		!温度
double precision::dt,dx								!刻み幅
double precision::Ndt,Ndx							!各軸の値を保存用
integer::i,j,k
double precision::r

dt=dble(time)/dble(Nt)
dx=dble(size)/dble(Nx)

r=TempCond*dt/dx/dx

allocate(U(0:Nx,0:Nt))

open(10,file="PPDEquation_FDM.d")
!-------------------初期値を設定----------------------
do i=1,Nx-1
	U(i,0)=0.0d0
end do
!-------------------境界条件を設定--------------------
do j=0,Nt
	U(0,j)=0.0d0
	U(Nx,j)=100.0d0
end do
!--------------------陽解法--------------------
	!-----------計算-------------------
do j=1,Nt
	do i=1,Nx-1
		U(i,j)=(1.0d0-2.0d0*r)*U(i,j-1)+r*(U(i+1,j-1)+U(i-1,j-1))
	end do
end do

	!-------------値を書きだす------------
do j=0,Nt
	Ndt=j*dt
	do i=0,Nx
		Ndx=i*dx
		write(10,*) Ndx,Ndt,U(i,j) 
	end do
	write(10,*)
end do

close(10)
stop
end program PPDE_FDM