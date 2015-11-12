!************************
!ガウス―ジョルダンの消去法
!20150930
!ver0.0
!************************
!AugMatrix::拡大係数行列(Augmented Matrix)
!
!例として
!x_1+x_2-2x_3=3
!5x_1+2x_2+x_3=1
!x_1-4x_2+3x_3=8
!を解いている
!

program Gauss_Jordan
implicit none
!------------------変数を定義----------
integer::N 								!何元の連立方程式を解くか,AugMatrixの行になる
integer::M 								!AugMatrixの列になる
double precision,dimension(:,:),allocatable::AugMatrix 			!拡大係数行列
double precision::A11,A22,A33			!対角成分、元の数だけ必要になる
double precision::Ai,Ak					!計算で必要な列成分
integer::i,j,k 							!doループ用変数
!-------------------------------------
N=3										!N元連立方程式に対応し、拡大係数行列の行数になる
M=N+1									!拡大係数行列の列数になる
allocate(AugMatrix(N,M))				!拡大係数行列の大きさを確定

!拡大係数行列の初期値
AugMatrix(1,1)=1.0d0
AugMatrix(1,2)=1.0d0
AugMatrix(1,3)=-2.0d0
AugMatrix(1,4)=3.0d0
AugMatrix(2,1)=5.0d0
AugMatrix(2,2)=2.0d0
AugMatrix(2,3)=1.0d0
AugMatrix(2,4)=1.0d0
AugMatrix(3,1)=1.0d0
AugMatrix(3,2)=-4.0d0
AugMatrix(3,3)=3.0d0
AugMatrix(3,4)=8.0d0

do k=1,N
	Ak=AugMatrix(k,k)
	do j=1,M
		AugMatrix(k,j)=AugMatrix(k,j)/Ak
	end do
	do i=1,N
		if(i.ne.k)then
		Ai=AugMatrix(i,k)
		do j=1,M
			AugMatrix(i,j)=AugMatrix(i,j)-AugMatrix(k,j)*Ai
		end do
		end if
	end do
end do


!A11=AugMatrix(1,1)
!do j=1,M
!	AugMatrix(1,j)=AugMatrix(1,j)/A11
!end do
!do i=1,N
!	Ai=AugMatrix(i,1)
!	if(i.ne.1) then
!		do j=1,M
!			AugMatrix(i,j)=AugMatrix(i,j)-AugMatrix(1,j)*Ai
!		end do
!	end if
!end do
!
!
!do j=2,M
!	AugMatrix(2,j)=AugMatrix(2,j)/A22
!end do
!do i=1,N
!	Ai=AugMatrix(i,2)
!	if(i.ne.2) then
!		do j=2,M
!			AugMatrix(i,j)=AugMatrix(i,j)-AugMatrix(2,j)*Ai
!		end do
!	end if
!end do
!
!A33=AugMatrix(3,3)
!do j=3,M
!	AugMatrix(3,j)=AugMatrix(3,j)/A33
!end do
!do i=1,N
!	Ai=AugMatrix(i,3)
!	if(i.ne.3) then
!		do j=3,M
!			AugMatrix(i,j)=AugMatrix(i,j)-AugMatrix(3,j)*Ai
!		end do
!	end if
!end do

!行列の中身を確認----------------
do i=1,N
	write(6,*) AugMatrix(i,1),AugMatrix(i,2),AugMatrix(i,3),AugMatrix(i,4)
end do


stop
end program Gauss_Jordan
