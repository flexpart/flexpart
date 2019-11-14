! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine cbl(wp,zp,ust,wst,h,rhoa,rhograd,sigmaw,dsigmawdz,tlw,ptot,Q,phi,ath,bth,ol,flagrein)
!                i  i i  i   i  i    i     i       i         i   o   o   o   o    o  i    o
     
    use par_mod, only:pi
    use com_mod, only:ldirect
    
    implicit none
!=======================================================================================================================================================
!=============== CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3 from LHB 2000                                      ========
!=============== LHH formulation has been modified to account for variable density profiles and backward in time or forward in time simulations ========
!=============== see Cassiani et al. BLM 2014 doi  for explanations and references                                                              ========
!=======================================================================================================================================================

    real :: usurad2,usurad2p,C0,costluar4,eps 
    parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=3,costluar4=0.66667,eps=0.000001)

       
    integer :: flagrein
    real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,rhoa,rhograd
    real ::fluarw,fluarw2
    real ::w3,w2
    real ::dw3,dw2
    real ::wb,wa
    real ::deltawa,deltawb
    real ::wold,wold2,wold_z
    real ::pa,pb,alfa
    real ::Phi,Q,ptot  
    real :: timedir
    real ::cuberoot
    real ::z0,ol,transition   
    

    real :: &
    erf, &
    aperfa, &
    aperfb, &
    ath, &
    bth

    real ::  &
    pow, &
    z, &    
    skew, &
    skew2, &
    radw2, &
    rluarw, &
    xluarw, &
    aluarw, &
    bluarw, &
    sigmawa, &
    sigmawb, &
    dskew, &
    dradw2, &
    dfluarw, &
    drluarw, &
    dxluarw, &
    daluarw, &
    dbluarw, &
    dsigmawa, &
    dsigmawb, &
    dwa, &
    dwb, &
    sigmawa2, &
    sigmawb2
    
    
    dens=rhoa      !put to 1 for test constant density simulation
    ddens=rhograd  !put to 0 for test constant density simulation
             
    
    timedir=ldirect !ldirect contains direction of time forward (1) or backward(-1)
    !========================= assegnazione z ===========================================================
    z=(zp/h)
    
    !================== stability transition function see Cassiani et al(2015) BLM ======================
    transition=1.
    !if (ol.lt.-50) transition=((sin(((ol+100.)/100.)*pi))-1.)/2.
    if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5

    !========================= momento secondo ==========================================================
    
    w2=(sigmaw*sigmaw)
    dw2=(2.*sigmaw*dsigmawdz)
    

    !=================== dissipazione =======================================
    
    alfa=2.*w2/(C0*tlw)

    !========================================================================
    
    wold=timedir*wp
    
    !=========================================================================
    !------------------------------ momento terzo ============================
                                
    w3=((1.2*z*((1.-z)**(3./2.)))+eps)*(wst**3)*transition
    dw3=(1.2*(((1.-z)**(3./2.))+z*1.5*((1.-z)**(1./2.))*(-1.)))*(wst**3)*(1./h)*transition
    
    !===========================================================================0
   
    skew=w3/(w2**1.5)
    skew2=skew*skew
    dskew=(dw3*w2**(1.5)-w3*1.5*w2**0.5*dw2)/w2**3
    radw2=w2**0.5
    dradw2=0.5*w2**(-0.5)*dw2
    !costluar4=0.66667  ! costante da LHH
    fluarw=costluar4*(cuberoot(skew))	!skew**(1./3.)
    fluarw2=fluarw*fluarw
    
    if (skew.ne.0) then
      
        dfluarw=costluar4*(1./3.)*cuberoot(skew**(-2.))*dskew
    
        rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
        xluarw=(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)    !----> r^1/2
        
        drluarw=( ((3.*(1.+fluarw2)**2*(2.*fluarw*dfluarw)*skew2)+ &
            (1.+fluarw2)**3*2.*skew*dskew) *(3.+fluarw2)**2.*fluarw2 - &
            (1.+fluarw2)**3*skew2* &
            ( (2.*(3.+fluarw2)*(2.*fluarw*dfluarw)*fluarw2) + &
            (3.+fluarw2)**2*2.*fluarw*dfluarw) )/ &
            (((3.+fluarw2)**2.*fluarw2)**2)
                                     
        dxluarw=( ((1.5*(1.+fluarw2)**0.5*(2.*fluarw*dfluarw)*skew)+ &
            (1.+fluarw2)**1.5*dskew) *(3.+fluarw2)*fluarw - &
            (1.+fluarw2)**1.5*skew* &
            (3.*dfluarw+3*fluarw2*dfluarw) )/ &
            (((3.+fluarw2)*fluarw)**2)
    
    else        
        dfluarw=0.
        rluarw=0.
        drluarw=0.
        xluarw=0.
        dxluarw=0.
    end if



   aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
   bluarw=1.-aluarw

   daluarw=-0.5*(  (dxluarw*(4.+rluarw)**0.5) - &
           (0.5*xluarw*(4.+rluarw)**(-0.5)*drluarw) ) &
           /(4.+rluarw) 
   dbluarw=-daluarw
   
   sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
   sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

   dsigmawa=dradw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5+ &
            radw2*( &
            (0.5*(bluarw/(aluarw*(1.+fluarw2)))**(-0.5))      * &
            ( &
            (dbluarw*(aluarw*(1.+fluarw2))- &
            bluarw*(daluarw*(1.+fluarw2)+aluarw*2.*fluarw*dfluarw)) &
            /((aluarw*(1.+fluarw2))**2) &
            ) &
            )
  dsigmawb=dradw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5+ &
           radw2*( &
           (0.5*(aluarw/(bluarw*(1.+fluarw2)))**(-0.5))      * &
           ( &
           (daluarw*(bluarw*(1.+fluarw2))- &
           aluarw*(dbluarw*(1.+fluarw2)+bluarw*2.*fluarw*dfluarw)) &
           /((bluarw*(1.+fluarw2))**2) &
           ) &
           )
    
           wa=(fluarw*sigmawa)
           wb=(fluarw*sigmawb)
           dwa=dfluarw*sigmawa+fluarw*dsigmawa
           dwb=dfluarw*sigmawb+fluarw*dsigmawb

           deltawa=wold-wa
           deltawb=wold+wb
           wold2=wold*wold
           sigmawa2=sigmawa*sigmawa
           sigmawb2=sigmawb*sigmawb
                            
           if (abs(deltawa).gt.6.*sigmawa.and.abs(deltawb).gt.6.*sigmawb) flagrein=1

           pa=(usurad2p*(1./sigmawa))*(exp(-(0.5*((deltawa/sigmawa)**2.))))
           pb=(usurad2p*(1./sigmawb))*(exp(-(0.5*((deltawb/sigmawb)**2.))))
                            	    
           ptot=dens*aluarw*pa+dens*bluarw*pb
                               
           aperfa=deltawa*usurad2/sigmawa
           aperfb=deltawb*usurad2/sigmawb

           Phi=-0.5* &
               (aluarw*dens*dwa+dens*wa*daluarw+aluarw*wa*ddens)*erf(aperfa) &
               +sigmawa*(aluarw*dens*dsigmawa*(wold2/sigmawa2+1.)+ &
               sigmawa*dens*daluarw+sigmawa*ddens*aluarw+ &
               aluarw*wold*dens/sigmawa2*(sigmawa*dwa-wa*dsigmawa))*pa &
               +0.5* &
               (bluarw*dens*dwb+wb*dens*dbluarw+wb*bluarw*ddens)*erf(aperfb) &
               +sigmawb*(bluarw*dens*dsigmawb*(wold2/sigmawb2+1.)+ &
               sigmawb*dens*dbluarw+sigmawb*ddens*bluarw+ &
               bluarw*wold*dens/sigmawb2*(-sigmawb*dwb+wb*dsigmawb))*pb

               Q=timedir*((aluarw*dens*deltawa/sigmawa2)*pa+ &
               (bluarw*dens*deltawb/sigmawb2)*pb)

               ath=(1./ptot)*(-(C0/2.)*alfa*Q+phi)
               bth=sqrt(C0*alfa)
              !bth=sngl(sigmaw*sqrt(2.*tlw))

    return
    end
                            




    FUNCTION CUBEROOT (X) RESULT (Y)

    IMPLICIT NONE

    real, INTENT(IN) :: X
    real:: Y

    real, PARAMETER :: THIRD = 0.333333333


    Y = SIGN((ABS(X))**THIRD, X)

    RETURN

    END FUNCTION CUBEROOT
    
    
    

    FUNCTION CUBEROOTD (X) RESULT (Y)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION :: Y

    DOUBLE PRECISION, PARAMETER :: THIRD = 0.33333333333333333333333333333333333333333333333333333333333333333333333333333333333D0


    Y = SIGN((ABS(X))**THIRD, X)

    RETURN

    END FUNCTION CUBEROOTD
