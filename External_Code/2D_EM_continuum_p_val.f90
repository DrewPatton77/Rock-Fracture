! Perform a Maximum Likelhood Estimation of Power-Law exponent by the
! false position method for all possible thresholds inside de data
!
!
! How to Use:
!
! willing to fit the exponent of a law : p(x) = x^t / N(t,E_low,E_high) 
!
! 1) build the file "whatever_name".Sorted containing the dataset in 
! a single column; sorted from smaller to bigger value  (you may consider 
! a quicksort algorithm to save time in long series)
!
! 2) Choose a "c" value (in this code) according to the density 
! of the grid you want to obtain. Make sure the size of the vector "f" is
! big enough for the data sample
!
! 3) execute the code:  2D-EM_continuum whatever_name  . The program will
! automatically extract information about the boundaries of the sample
!
! 4) the output data is stored in files "whatever_name".Ctr and "whatever_name".EMap.
! the column format is:
! E_low     E_max     fitted exponent     error from lik.function    n of values inside this interval 
! "whatever_name".Ctr  is intended to help other programs to calculate 
! the contour lines for the map.
! "whatever_name".Emap makes easy to represent the E_low cuts for diferent
! E_high values, similarly to Clauset's analysis.
! Both of them work well for representing the maps
! "whatever_name".Histo stores the frequency of each exponent. Just recomended for
! quick checks before map representation.
! "whatever_name".CCDF representents the Complementary Comulative Distributioon
! Function for the amusement of the reader.
!
! (5) p-values or other statistical tests are highly recomended in order 
! to test the validity of the results
!
! author: Jordi Baró i Urbea, EMC dept., Universitat de Barcelona
! ref: Phys. Rev. E 85, 066121
! http://pre.aps.org/abstract/PRE/v85/i6/e066121

program ExponentMaps
implicit real*8(a-z)
real (kind=8), parameter :: tau0 = 1.0D0 
real (kind=8), parameter :: tau1 = 3.5D0

! per fer els histogrames
integer(kind=4), parameter :: tauN = 80
real(kind=8), parameter    :: tauBin = (tau1 - tau0)/tauN
real (kind = 8) :: tauHisto(tauN), tauHistoLin(tauN),tauTrue, bin, NhT, NhTl

!integer (kind=4)    :: C,D
integer (kind=8) :: i,j,Amax_max
integer (kind=4) :: istat=0
real (kind=8) ::  Amin, Amax, Amax2, Amin2

real(kind=8) ::  A_thres
real (kind = 8) :: tauMean, tauDes, sum, nAllau
integer(kind=8) nTaus
!integer, pointer :: pf(:)

!#######################################################################
!+++++++++++++++++++ define threshold step +++++++++++++++++++++++++++++
real (kind=4) :: c=1.1 ! next step x(i+1) = c*x(i)
!+++++++++++++++++++ vector of data f with redifinable max. dimension 
real*8 f(17000000) 
!#######################################################################

common f, tauTrue, Amax_max
!real*8 t,P,Q,x
!real*8 DL


character*64 fname

!R='2.17'
!L='256'
CALL GETARG(1 , fname)
!CALL GETARG(2 , R)

print *, '=============================================================='
print *, '2D-PowerLaw Exponent Mapped in Amax and Amin for Continum data'
print *, '=============================================================='
print *, ''

print *, 'reading data...'
!llegeix el fitxer i assigna i i f[i]
open (unit=11,file=trim(fname)//'.Sorted', status = 'old')
open (unit=13,file=trim(fname)//'.EMap', status = 'unknown')
open (unit=14,file=trim(fname)//'.Histo', status = 'unknown')
open (unit=15,file=trim(fname)//'.CCDF', status = 'unknown')
open (unit=16,file=trim(fname)//'.Ctr', status = 'unknown')
open (unit=17,file=trim(fname)//'.Xmin', status = 'unknown')
!open (unit=10, file='taus/tottaus19', status='unknown', position = 'append')
!L = 32
!Amax_max = 64*64*64

nAllau=0
i=1;
do while (istat.EQ.0);
 read (11,*,iostat=istat) f(i)
 nAllau = nAllau + f(i)
 i=i+1;
! print*,f(i)
end do
Amax_max = i-2
close (unit=11)


print *, 'prepare comulative (optional)'
sum = 0
do i=Amax_max, 1, -1
  sum = sum + 1.0D0/Amax_max
  write(15,*) f(i), sum
end do
close (unit = 15)


 print *, 'procesing exponents...(', f(1) , f(Amax_max) , ')'

! versió especial
!A_thres = 1.0D8
!normalment
A_thres = c*f(Amax_max)
!
Amax=f(1)*5
do while (Amax2.lt.A_thres)
Amin=f(1)
! Amax2 = ceiling(Amax*1.1)
 Amax2 = Amax*c

 do while(Amin.lt.Amax .and. Amin .lt.f(Amax_max))  
   Amin2 = Amin*c

       call tau( Amin, Amax,.false.)

    !	promitjos   

    !	histograma   
       bin = tau0
       j = 0
       if ((tauTrue .le. tau1) .and. (tauTrue.ge.tau0) ) then
        do while (j.lt.tauN .and. tauTrue.gt.bin)
           j = j+1
           bin = bin + tauBin
        end do
        tauHisto(j)= tauHisto(j)+1
        tauHistoLin(j) =tauHistoLin(j) + (Amin2-Amin)*(Amax2-Amax)
        nTaus = nTaus + 1
        tauMean = tauMean + tauTrue
        tauDes = tauDes + tauTrue*tauTrue
       endif
    Amin=Amin2
  !      print*, 'aMax = ', Amax,'aMin = ', Amin, ' tau = ', tauTrue
   end do

  Amax=Amax2
  write(13,*) ''
  write(13,*) ''
!  print*, nTaus, tauMean, tauDes
 end do
! close(unit = 13)

Amax= A_thres
Amin=f(1)

do while(Amin.lt.Amax .and. Amin .lt.f(Amax_max))  
   call tau( Amin, Amax, .true.)
   Amin = Amin*c
end do

print*,'processant distribution of exponents (optional)...'
NhT = 0
NhTl = 0
do i=1,tauN
  NhT = NhT + tauHisto(i)
  NhTl = NhTl + tauHistoLin(i)
!  print*, NhT
end do
bin= tau0+0.5*tauBin
do i=1,tauN
  bin = bin + tauBin
  write(14,*) bin, tauHisto(i)/NhT, tauHistoLin(i)/NhTl
end do




!tauMean = tauMean/nTaus
!tauDes = (tauDes/nTaus - tauMean*tauMean)
!tauTrue = sqrt(tauDes/(nTaus-1))

!write(10,*) L//'  '//R, tauMean, tauDes, tauTrue


close(unit = 13)
close(unit = 14)
close(unit = 15)
close(unit = 16)
close(unit = 17)
end program


subroutine tau( Amin, Amax, last)
! Calcula la tau per tots els amin amax


 implicit real*8(a-z)
 
 logical last
 real*8 Amin, Amax,  lmin, lmax, A_thres
 real*8 tauTrue, DL1,DL2,DL0,t0,t1,t2
!real (kind=8), parameter :: tau0,tau1
 real*8 f(17000000)
  real*8 epsilon
real*8 p_val, kmax

integer (kind=8) :: n, i, i0, Amax_max,j
 real*8 Q,P
common f, tauTrue, Amax_max, A_thres
! real*8 P
! P=0.0D0


  epsilon = 5.0D-3
t0=-1.9
t1=1.3

i=1

 do while (f(i) < Amin) 
   i=i+1
 end do
 Q=0.0D0
 i0 = i
 do while (f(i) < Amax .and. i.le. Amax_max)
    Q = Q + log(f(i))
    i=i+1
 end do

 n = i-i0
 
 if (n>0) then
 !Q=Q/n
lmin = log(Amin)  
lmax = log(Amax)
! print*, lmax
 DL0=DL(t0, Amin, Amax, lmin ,lmax,Q,n)
 DL1=DL(t1, Amin, Amax, lmin ,lmax,Q,n)
! if (DL1*DL0 > 0) then
!   print*,'cagada', DL1 , DL0
!   tautrue = 0
!   return
! endif
 
 DL2=5
! print*, DL0, DL1
 do while (abs(t1-t0)>epsilon)
!print*, 'Its a trapp!!', DL0, DL1, t2
!calcula el L' en base a la nove tau
    t2=t1 - DL1*(t1-t0)/(DL1-DL0)
    DL2 = DL(t2, Amin, Amax,lmin ,lmax,Q,N)
 
 ! if (DL1*DL2<0) then 
 !  t0  = t2 
 !  DL0 = DL2
 ! else
 !  t1  = t2 
 !  DL1 = DL2
 ! end if
 
     t0=t1
     t1=t2
     DL0=DL1
     DL1=DL2
 end do
  
! print*,Q, n, Amin, Amax
 
 tauTrue = 1.-t2
 
! print*,'here', i, i0, tauTrue, t2

! calcula p-vals
 kmax = 0
 P = (Amax/Amin)**(t2)
 Q= 1.0D0/(1.0D0-P)


! print*,'here', tauTrue, t2, Q,P

 do j=i0, i-1
     kmax = max(kmax,  abs(Q*((f(j)/Amin)**(t2) - P) -1 + (j-i0)/(real(n))))
 end do
 p_value = P_val((kmax),sqrt(real(n)),150)
 
 if(last) then
   write(17,*)  Amin,Amax, tauTrue, sqrt(abs((t0-t1)/(DL0-DL1))), n, kmax, p_value
 else
   write(13,*)  Amin,Amax, tauTrue, sqrt(abs((t0-t1)/(DL0-DL1))), n, kmax, p_value 
   write(16,*)  Amin,Amax, tauTrue, sqrt(abs((t0-t1)/(DL0-DL1))), n, kmax, p_value
 end if
end if 

 
 return 
end subroutine tau


real*8 function DL(t,Amin, Amax, lmin, lmax, Q, N)
implicit real*8(a-z)
! calcula el DL per falsa posició
 real*8 lmin,lmax, q,t, Amin,Amax
 integer*8 N
! integer f(100000)
 DL = N*(((Amax**(t))*lmax -(Amin**(t))*lmin)/((Amax**(t)) - (Amin**(t))) - 1./(t)) - Q

 return
end 


!real*8 function P_val(x, kmax)
real*8 function P_val(x, sqrn, kmax)
! prob(K<=x) sqrt(2*pi)/x * sum_k (exp(-(2*k-1)*(2*k-1)*pi*pi/(8*x*x) )
!sqrt(2*pi) = 2.506628275
!pi*pi/8 = 1.23370055

   real (kind=8) :: x,  D
   real (kind=4) :: sqrn
   integer (kind=4) :: k, kmax

   D=0
   do k=1, kmax
!    D= D+ exp(-(2*k-1)*(2*k-1)*1.23370055/(x*x) )
!    D= D+ exp(-(2.0*k-1)*(2.0*k-1)*1.23370055/(x*x) )
!      D= D - (1.0D0-2.0D0*mod(k,2))*exp(-2.0D0*k*k*x*x)
      D= D + (- 1.0D0)**(k-1) *exp(-2.0D0*(k*(x*sqrn+0.12*x+0.11*x/sqrn))**2)
   end do 
!   P_val = 1.0D0-2.506628275*D/x
   P_val = 2.0D0*D
 
   return 
end
