program SVD

implicit none


integer :: i,j,k,lx,m,n,Ns,Nx,fldr,Nt,Nft,lz,ls,lt
real, allocatable :: Daux(:,:),DMNO(:,:),D1(:,:),Xaux(:),Daaux(:,:),A(:),env(:),Aaux(:,:)
integer, allocatable:: count2(:),count1(:)

Ns=1751 ! Numero de amostras
Nx=98514 ! numero de traços
Nt=117 !Numero de receptores
lx=3 !Janela Movel em x
lz=3 !Janela movel em z
fldr=841 !numero de tiros

ls=ns+2*lz
lt=nt+2*lx

allocate (Daux(lz,lx),D1(Ns,Nt),count2(lt),DMNO(Ns,Nt),Xaux(ls),Daaux(Ns,Nt),A(ls), count1(ls))

D1=0.
Daux=0.
DMNO=0.


	call system ('rm copex.ad')

       open(10,file='l2140266-geom-ed2-water.ad',form='unformatted',status='unknown',recl=4*Ns,access='direct')
       open(20,file='copex.ad',form='unformatted',status='unknown',recl=4*Ns,access='direct')
       open(40,file='N.txt',status='unknown')

       
do m=1,fldr

        daaux=0.
        DMNO=0.

         do k=1,nt
                write(*,*) "tiro", m,'receptor',k

                read(10,rec=Nt*(m-1)+k)Daaux(1:ns,k)
        end do
     
        call LMNO(Ns,Nt,Daaux(1:ns,:),DMNO(1:ns,:))
            
        count2=0.;count1=0.;D1=0.
 

 call svd_filt(Ns,nt,lz,lx,DMNO,D1)

               Do i=1,nt
                      ! write(20,rec=Nt*(m-1)+i)DMNO(:,i)
                      write(20,rec=Nt*(m-1)+i)DMNO(:,i)-D1(:,i)

               end do

end do


Deallocate (Daux,D1,count2,DMNO,Daaux,Xaux)


call system ('xwigb n1=1751 perc=90  title=SVD legend=1 <copex.ad &')

 close(10)
 close(20)

stop

end program SVD
	

 !##############################################################
        subroutine svd_filt(nz,nx,lz,lx,D,Dfilt)
        implicit none
        integer         :: i,j,nz,nx,lz,lx
        real            :: D(nz,nx),Daux(1-lz/2:nz+lz/2,1-lx/2:nx+lx/2),Dfilt(nz,nx),D1(lz,lx)
                
        Daux = 0.
        Daux(1:nz,1:nx) = D

        do j=1,nx
                do i=1,nz
                        call svd_power(lz,lx,Daux(i-lz/2:i+lz/2,j-lx/2:j+lx/2),D1)
                        Dfilt(i,j)=D1(lz/2+1,lx/2+1)

                end do        
        end do

        return

        end subroutine       


!######################################################################
!============================================================ 
        subroutine svd_power(m,n,X,X1)
!============================================================
! devolve em X1 a primeira autoimagem de X
!Subrotina para utilizar o power method SVD
! use n impar
!09/Fevereiro/2013
!Rafael Rodrigues Manenti
!m e n são as ordens da matriz
!A é a matriz de entrada
!Lambda é o autovalor dominante ou primeiro autovalor
!x é o autovetor dominante ou primeiro autovetor
!n_iter é o número de iterações utilizadas para se achar o valor aproximado de lambda e x.
!n_iter está fixado em 50 iterações
!-----------------------------------------------------------------------
! modifificada Milton J. Porsani 12/02/2013 (terça de carnaval)
! modificada MJP p/ devolver apenas a primeira autoimagem 05/09/2018
!============================================================
      integer :: m,n,ikey,niter_max,k
      real :: X,X1,XTX,W1,Xaux,U,V,tol,Q0,Q1,Sigma,XX,Slambda


      dimension  X(m,n),X1(m,n)
      allocatable XTX(:,:),w1(:),xaux(:),u(:),v(:)
      niter_max=50 ; tol=1.e-30     
      allocate (XTX(n,n),w1(n),xaux(n),u(m),v(n))
!------------------------------------------------------
      u=0.;v=0.
      if(sum(x*x).eq.0.)return
      xtx=matmul(transpose(X),X)
      xaux=1.
      w1=matmul(xtx,xaux)
      Q0=dot_product(xaux,w1)
      xaux=w1/sqrt(Q0)
!========================================================
     ikey=0
     k=0
     do while(ikey.eq.0)
        k=k+1                          
         w1=matmul(xtx,xaux)
         Q1=dot_product(xaux,w1)      
         xaux=w1/sqrt(Q1)                     
         if(abs(Q1-Q0).lt.tol.or.k.eq.niter_max)ikey=1  
         Q0=Q1                              
      end do                                     
!=========================================================        
      xx=dot_product(xaux,xaux)                                                            
      slambda=Q1/xx              ! sp positivo                                             
      sigma=sqrt(slambda)       !
!=========================================================
      V(1:n)=xaux(1:n)/sqrt(xx)      !normalizando p/ ter energia unitária                  
      u=matmul(x,v) 
      
      do k=1,m
        x1(k,1:n)=v(1:n)*u(k)
      end do
!========================================================
      deallocate (XTX,w1,xaux,u,v)
      
      return
      end subroutine

subroutine LMNO(Ns,Nt,A,DMNO)


integer :: Ns,Nt,V
real :: T0,dx,dt
real, Allocatable ::T(:)
real, dimension(Ns,Nt) ::DMNO,A
allocate (t(Ns))
nt=117



Xaux=0


V=1550
t0=0
dx=25!distancia entre receptores
dt=0.004!taxa de amostragem

do i=Nt,1,-1

	t(i)=(t0 +(((i-1)*dx)/v))

	dNt=t(i)/dt

	
	N=NINT(dNt)

	write(40,*)N
        
	do k=1,Ns-N

	DMNO(k,Nt-i+1)=A(k+N,Nt-i+1)

	end do
	



	
end do

deallocate (t)



return 

end





