!**********************
!2次元ラプラス方程式を差分法で
!Gauss-Saidel法と逐次緩和法(Successive over-relaxation method:SOR法)
!を用いて解く
!20151002
!ver0.0
!**********************

program SORM
implicit none

!------------------変数を定義---------------------------
double precision,parameter::Omega=1.5d0 				!Omega=1でガウスザイデル法
integer,parameter::m=1									!x軸の境界0-m
integer,parameter::n=1									!y軸の境界
integer,parameter::h=10									!x方向の刻み数
integer,parameter::l=10									!y方向刻み数
double precision,parameter::Converge=0.0001d0			!収束条件
integer,parameter::Max=100								!収束するまでの最大計算回数
double precision,dimension(:,:),allocatable::U			!U(x,y)
double precision,dimension(:,:),allocatable::LastU		!SOR法収束判定用前のU(x,y)
integer::i,j,k											!doループ用
double precision::dm,dn 								!刻み幅
double precision::U1
double precision::Sum0,Sum
double precision::x,y
allocate(U(0:h,0:l))
allocate(LastU(0:h,0:l))

dm=dble(m)/dble(h)		!x刻み
dn=dble(n)/dble(l)		!y刻み

open(10,file="SORM.d")

!------------------境界条件----------------------------
do j=0,l
	U(0,j)=0.0d0
	U(h,j)=0.0d0
end do
do i=0,h
	U(i,0)=0.0d0
	U(i,l)=4*i*dm*(1-i*dm)
end do
!-------------------初期値を設定(境界条件以外)------------
do i=1,h-1
	do j=1,l-1
		U(i,j)=0.0001d0
	end do
end do
!-----------------------SOR法---------------------------
do k=1,Max
	Sum=0.0d0				!初期化
	Sum0=0.0d0				!初期化

	do i=1,h-1
		do j=1,l-1
			LastU(i,j)=U(i,j)
			U1=(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))/4.0d0-U(i,j)
			U(i,j)=U(i,j)+Omega*U1

			Sum0=Sum0+abs(U(i,j))				!収束判定用に和をとっておく
			Sum=Sum+abs(U(i,j)-LastU(i,j))		!収束判定用の計算
		end do
	end do
	
	!-----------収束判定-----------
	Sum=Sum/Sum0
	if(Sum<=Converge) exit

end do

do i=0,h
	do j=0,l

		x=i*dm
		y=j*dn

		write(10,*) x,y,U(i,j)
	end do
end do



deallocate(U,LastU)		!allocateで確保してたのを解放
close(10)
stop
end program SORM