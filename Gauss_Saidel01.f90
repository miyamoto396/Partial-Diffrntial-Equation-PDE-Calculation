!*********************
!ガウス-ザイデル(Gauss-Saidel)の反復法で連立方程式をとく
!20151001
!ver0.1 
!収束条件を相対誤差にかえてみる
!結果収束しなくなった。だめじゃん
!*********************
!Coeff(N:N):係数行列Coefficien
!Const(N):定数項配列constant
!X(N)
!LastX(N):収束判定用の配列。前のX(N)を格納する

program Gauss_Saidel
implicit none
!------------------変数を定義-----------------
integer,parameter::N=3									!N元連立方程式に対応する
integer,parameter::Max=100								!収束するまでの最大計算回数
double precision,parameter::Converge=0.00001d0			!収束条件
double precision,dimension(:,:),allocatable::Coeff 		!係数行列
double precision,dimension(:),allocatable::Const		!定数項配列
double precision,dimension(:),allocatable::X			!変数xの配列
double precision,dimension(:),allocatable::LastX		!収束判定用
double precision::Judge									!収束判定用に用意した変数=0なら収束
integer::i,j,k
double precision::Sum
!-------------------------------------------

allocate(Coeff(N,N))
allocate(Const(N))
allocate(X(N))
allocate(LastX(N))

!係数行列の値
Coeff(1,1)=2.0d0
Coeff(1,2)=1.0d0
Coeff(1,3)=1.0d0
Coeff(2,1)=2.0d0
Coeff(2,2)=3.0d0
Coeff(2,3)=1.0d0
Coeff(3,1)=1.0d0
Coeff(3,2)=1.0d0
Coeff(3,3)=3.0d0

!定数項の値
Const(1)=2.0d0
Const(2)=4.0d0
Const(3)=-1.0d0

!Xの初期値の設定
do i=1,N
	X(i)=0.0d0
end do

!---------------------反復法の計算-----------------
do k=1,Max
	Judge=0.0d0			!初期化
	Sum=0.0d0
	do i=1,N 	   		!収束判定用の一個前の解を入れておく
	LastX(i)=X(i)
	end do

	do i=1,N 			!X(i)を求める
		do j=1,N
			if(i.ne.j)then
				Sum=Sum+Coeff(i,j)*X(j)
			end if
		end do

		X(i)=(Const(i)-Sum)/Coeff(i,i)
	end do

	!---------収束判定--------------
	do i=1,N
		if( abs((X(i)-LastX(i)/X(i))) > Converge)then
			Judge=Judge+1.0d0
		end if
	end do
	!------------------------------
	!write(6,*) Judge
	if(Judge.ne.0)cycle			!Judge/=0のときまた反復して計算する
	exit						!Judge=0なら計算終了
end do
!-------------------------------------------------

write(6,*) "最大計算回数Max=",Max
write(6,*) "実行計算回数k=",k					!実際に反復した回数
write(6,*) "収束条件Converge=",Converge		!収束条件

write(6,*) "-----------solution-----------"
do i=1,N
	write(6,*) X(i)
end do

end program Gauss_Saidel