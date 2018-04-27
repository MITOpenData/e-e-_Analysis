C***********************************************************************
C
C This is a set of routines to perform a standard QCD analysis of 
C event shapes and inclusive charged particle distributions
C AT LEP1 (!!) with ALPHA.
C Histogram-based, not Ntuple-based
C
C Note : for an analysis at LEP2 the event selection has to be
C        adopted correspondingly
C
C Author : Guenther Dissertori / 23-MAY-01
C
C Version: 1.0 First setup
C***********************************************************************

      SUBROUTINE Quinit
C-----------------------------------------------------------------------
C ALPHA user initialization routine
C =================================
C
C Author: Guenther Dissertori  
C-----------------------------------------------------------------------

      Implicit NONE 
#include "qdecl.h"
#include "qcde.h"

      Integer   i,j,Istat
 
C -- a counter of good hadronic events
      Integer   Nhad
      Common /Norm/ Nhad

C -- definition of histos
      Include 'histonaming.h'

#include "qmacro.h"
C-----------------------------------------------------------------------

C--      standard initialization

      call QCD_Init

C--      specific initalization

      Nhad = 0

C -- book histos
C -- note: histos with identifier below 40 are used for
C -- standard histos filled in the initialization!

C --NOTE : BINNING MIGHT BE OPTIMIZED FOR SPECIFIC NEEDS

C -- event shapes
      call HBOOK1(idth,'Thrust       ',50,0.5,1.0,0.)
      call HBOOK1(idta,'Thrust Major ',50,0.0,1.0,0.)
      call HBOOK1(idtm,'Thrust Minor ',50,0.0,0.5,0.)
      call HBOOK1(idmh,'Mh2/s        ',50,0.0,0.5,0.) 
      call HBOOK1(idcp,'Cparameter   ',50,0.0,1.0,0.)
      call HBOOK1(idbt,'Btot         ',50,0.0,0.4,0.) 
      call HBOOK1(idbw,'Bwide        ',50,0.0,0.4,0.) 
      call HBOOK1(idly,'-lny3        ',50,0.0,10.,0.) 

C -- D parameter for different topologies
      call HBOOK1(iddi,'Dpar, incl   ',50,0.,1.,0.)   
      call HBOOK1(idd1,'Dpar, y=0.1  ',50,0.,1.,0.)   
      call HBOOK1(idd2,'Dpar, y=0.05 ',50,0.,1.,0.)   
      call HBOOK1(idd3,'Dpar, y=0.025',50,0.,1.,0.)   

C -- moments
      call HBOOK1(idmo,'Moments of Thrust',8,0.5,8.5,0.)

C -- three-jet rate for three ycut values
      call HBOOK1(idjr,'3 jets',3,0.5,3.5,0.)

C -- angular correlations in 4-jet events
      call HBOOK1(idbz,'|cos(BZ)| ',20,0.,1.,0.)   
      call HBOOK1(idnr,'|cos(NR)| ',20,0.,1.,0.)   
      call HBOOK1(idks,'cos(KSW)  ',20,-1.,1.,0.)   
      call HBOOK1(id34,'cos(a34)  ',20,-1.,1.,0.)   

C -- inclusive charged particle distributions
      call HBOOK1(idxp,'xp   ',50,0.,1.,0.)  ! = pch/pbeam
      call HBOOK1(idxi,'Xsi  ',50,0.,6.,0.)  ! = -ln(xp)
      call HBOOK1(idra,'y    ',50,0.,6.,0.)  ! rapidity w.r.t. thrust axis 
      call HBOOK1(idpi,'PTin ',50,0.,15.,0.) ! momentum in event plane 
      call HBOOK1(idpo,'PTout',50,0.,5.,0.)  ! momentum out of event plane

C -- angle of thrust axis
      call HBOOK1(idt1,'|cos(theta_thrust)| bef',100,0.,1.,0.) ! before selection
      call HBOOK1(idt2,'|cos(theta_thrust)| aft',100,0.,1.,0.) ! after selection
      call HBOOK1(idt3,'|cos(theta_thrust)| chg',100,0.,1.,0.) ! after selection,
C                                                                chg tracks only

C -- some info and couting
      call HBOOK1(idc1,'Iret',21,-20.5,0.5,0.)
      call HBOOK1(idc2,'counting',6,0.5,6.5,0.)

C -- some monitoring
      call HBOOK1(idet,'Etot sel.obj',150,0.0,150.0,0.)
      call HBOOK1(idnt,'Numb sel obj',100,0.0,100.0,0.)
      call HBOOK1(idfl,'flavour',4,-0.5,3.5,0.)


      CALL HMINIM (0,0.0)
      CALL HIDOPT (0,'STAT')
      CALL HIDOPT (0,'BLAC')
      CALL HIDOPT (0,'NPHI')
      CALL HSQUEZ ('YES')


      RETURN
      END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------



      SUBROUTINE Quevnt(qt,kt,qv,kv)
C-----------------------------------------------------------------------
C ALPHA Event Analysis Routine
C ============================
C
C Author: Guenther Dissertori 
C-----------------------------------------------------------------------

      Implicit NONE
#include "qdecl.h"
#include "qcde.h"
 
      Real       qt,qv
      Integer    kt,kv
C-----------------------------------------------------------------------

      Logical   passed,PassedSel
      Integer   iret,prob_code(10)
      Integer   ntrk,itkl,i,j,ITHRU,IretGlobal
      Integer   Scheme,Error,Nch,Ntr
      Integer   IFLAV,ISTOG,NLINK,JKEVH,Flavour
      External  NLINK

      Real      THRUST,Ech,Cts,Etot,px,py,pz


      Integer   Nhad
      Common /Norm/ Nhad

C -- histo definitions
      Include 'histonaming.h'

#include "qmacro.h"
C-----------------------------------------------------------------------


C -- very first initialization
      IretGlobal = 0
      passed     = .False.
      PassedSel  = .False.

C -- do counting of all events which are looked at
      call HFILL(idc2,1.,0.,1.)

C -- check ALEPH status
      call ALEPH_STATUS(iret,prob_code)
      if (iret.ne.0) then
        if (KNEVT.lt.1000) then
          write(6,'(A,I5,9I2)')
     &         'ALEPH_STATUS failed : KNEVT,prob_code = ',KNEVT,
     &         (prob_code(i),i=1,9)
        endif
        IretGlobal = -1
        call HF1(idc1,float(IretGlobal),1.)
        RETURN                  ! --------->     
      endif

C -- compute thrust axis with the preselected objects
      call QJOPTR('EF',' ')
      call QJTHRU(THRUST,'THRUST',KRECO)
      if (THRUST.eq.0.) then
        IretGlobal = -2
        call HFILL(idc1,float(IretGlobal),0.,1.)
        write(6,*) 'No thrust',KRUN,KEVT
        RETURN                ! ------------>
      endif
      ITHRU = KPDIR('THRUST',KRECO)
      if (ITHRU.le.0.or.QP(ITHRU).eq.0.) THEN
        CALL QWMESE('Problems in calculating thrust axis')
        IretGlobal = -3
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                ! ------------>
      endif
C -- store the angle of the axis before selection
      call HFILL(idt1,abs(QCT(ITHRU)),0.,1.)


C========================================================================
C -- select hadronic events according to standard procedure
C -- selection is done on detectorlevel!!!!
C -- PassedSel means w.r.t. to the simple selection criteria
C========================================================================
      call QCD_HadSel(KRECO,Nch,Ech,Cts,Ntr,PassedSel)
      if (.NOT.PassedSel) then
        IretGlobal = -4
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                  ! --------->
      endif

C -- we passed the selection, then go.....
C -- increase counter of hadronic events
      Nhad = Nhad+1
      call HFILL(idc2,2.,0.,1.)
      if (mod(Nhad,100000).eq.0.) then
        write(6,*)
        write(6,'(A,I8,A)')
     &       ' -----> ',Nhad,' events have been processed...'
      endif

C -- here additional event selection criteria might be put
      passed  = .True.

C========================================================================
C -- Select Enflw-Objects (additional track selection
C -- for calculation of event shapes)
C========================================================================
      call QCD_Enflow(KRECO,'ALL',ntrk,etot,px,py,pz,passed)

C -- monitoring
      call HFILL(idet,etot,0.,1.)
      call HFILL(idnt,float(ntrk),0.,1.)
      
C -- get the flavour
      Flavour = 0
C -- first, if MC
      if (XMCEV) then
        JKEVH=NLINK('KEVH',0)
        IF(JKEVH.EQ.0) THEN
          IFLAV=0
          ISTOG=0
        ELSE
          ISTOG=ITABL(JKEVH,1,4)
          IFLAV=ISTOG/1000
          IF(IFLAV.GE.10) IFLAV=IFLAV-10
        ENDIF
        if (IFLAV.eq.0) then
          CALL QWMESE('Event with no Flavour info')
        endif
        if     (IFLAV.lt.4) then
          Flavour = 1
        elseif (IFLAV.eq.4) then
          Flavour = 2
        elseif (IFLAV.eq.5) then
          Flavour = 3
        else
          CALL QWMESE('Do not know this flavour....')
        endif
      else
C -- else data : we don't do a flavour tag here, but could be added
      endif
      call HFILL(idfl,float(Flavour),0.,1.)

C -- compute again thrust axis with the selected objects and events
      call QJOPTR('EF',' ')
      call QJTHRU(THRUST,'THRUSTSEL',KRECO)
      if (THRUST.eq.0.) then
        IretGlobal = -5
        call HFILL(idc1,float(IretGlobal),0.,1.)
        write(6,*) 'After selection : No thrust',KRUN,KEVT
        RETURN                ! ------------>
      endif
      ITHRU = KPDIR('THRUSTSEL',KRECO)
      if (ITHRU.le.0.or.QP(ITHRU).eq.0.) THEN
        CALL QWMESE('After sel:Problems in calculating thrust axis')
        IretGlobal = -6
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                ! ------------>
      endif
C -- store the angle of the axis after selection
      call HFILL(idt2,abs(QCT(ITHRU)),0.,1.)


C -- Dump an event
      if (KNEVT.eq.2) then
        write(6,*)
        call QCD_EvDump
        write(6,*)
      endif

C========================================================================
C -- calculations of event shapes  + jet rates
C========================================================================
      call Calc_Evshapes(iret)
      
C -- bookkeeping of events for which all calculations ok
      if (iret.eq.0) call HFILL(idc2,5.,0.,1.)

C========================================================================
C -- calculations of inclusive charged particle distributions
C========================================================================
      call Calc_Inclusives(iret)

C -- bookkeeping of events for which all incl. calculations ok
      if (iret.eq.0) call HFILL(idc2,6.,0.,1.)

C -- fill status code
      IretGlobal = 0
      call HFILL(idc1,float(IretGlobal),0.,1.)

C -- that's it....

 9999 RETURN 
      END


C----------------------------------------------------------------------
C----------------------------------------------------------------------


      SUBROUTINE Quterm
C-----------------------------------------------------------------------
C ALPHA user termination routine
C =================================
C
C Author: Guenther Dissertori  / 5.1.95  
C-----------------------------------------------------------------------

      Implicit NONE
#include "qdecl.h"
#include "qcde.h"

      Integer      i,j

      Real         Ngood,N3j,N4j,N3e,N4e
      Real         sum1,sum2,mean,erro

      Integer   Nhad
      Common /Norm/ Nhad

      Real         Hi,Hie,Hstati,Hij
      Logical      Hexist
      External     Hi,Hie,Hexist,Hstati,Hij

C -- histo definition
      Include 'histonaming.h'

#include "qmacro.h"
C-----------------------------------------------------------------------


      write(6,*)
      write(6,*)'*****************************************************'
      write(6,*)' ====================> final calculations ....'
      write(6,*)


      if (Nhad.eq.0) then
        write(6,*)'No hadronic event on detectorlevel found !!!!!'
        RETURN
      endif

C -- number of accepted hadronic events
      write(6,*)
      write(6,'(1X,A,I8)')'Nhad used in normalization = ',Nhad
      write(6,*)
      
      write(6,'(2X,A,F10.0,A,F10.0,A)')'Out of ',HI(idc2,1),' events ',
     >      HI(idc2,2),' passed the cuts.'
      write(6,*)

      write(6,'(2X,A,F10.0,A,F10.0,A)')'Out of ',HI(idc2,2),' events ',
     >      HI(idc2,3),' had good event shape calculations.'
      write(6,*)

      write(6,'(2X,A,F10.0,A,F10.0,A)')'Out of ',HI(idc2,2),' events ',
     >      HI(idc2,4),' had good 4-jet calculations.'
      write(6,*)

      write(6,'(2X,A,F10.0,A,F10.0,A)')'Out of ',HI(idc2,2),' events ',
     >      HI(idc2,5),' had good event shape + jet calculations.'
      write(6,*)

      write(6,'(2X,A,F10.0,A,F10.0,A)')'Out of ',HI(idc2,2),' events ',
     >      HI(idc2,6),' had good inclusive chg tracks calculations.'
      write(6,*)
      
C -- compute the 3-jet rates
      Ngood = HI(idc2,2)
      if (Ngood.gt.0) then 
        do i=1,3
C -- get the numbers
          N3j = HI(idjr,i) 
          N3j = N3j/Ngood
          N3e = sqrt(N3j*(1-N3j)/Ngood)

          write(6,990)' 3-jet rate for ycut ',Ycut3j(i),' = ',N3j,N3e
        enddo
      endif


C -- compute the 4-jet rate
      Ngood = HI(idc2,2)
      if (Ngood.gt.0) then 
C -- get the numbers
        N4j = HI(idc2,4) 
        N4j = N4j/Ngood
        N4e = sqrt(N4j*(1-N4j)/Ngood)
        write(6,*)
        write(6,990)' 4-jet rate for ycut ',Ycut4j,' = ',N4j,N4e
      endif

C -- compute the first moment of thrust
      Ngood = HI(idc2,3)
      if (Ngood.gt.1) then 
C -- get the numbers
        sum1 = HI(idmo,1)
        sum2 = HI(idmo,5)
        mean = sum1/Ngood
        erro = Ngood/(Ngood-1.) * (sum2/Ngood - mean**2)
        erro = sqrt(max(1.e-7,erro/Ngood))
        write(6,*)
        write(6,991)' < 1 - Thrust > = ',mean,erro
      endif

      write(6,*)
      write(6,*)' ====================> end of final calculations ....'
      write(6,*)'*****************************************************'
      write(6,*)

 990  Format(2X,A,F6.3,A,F7.5," +- ",F7.5)
 991  Format(2X,A,F7.5," +- ",F7.5)

      RETURN
      END
      

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE ErrMess(Message,Var)
C-----------------------------------------------------------------------
C Routine which prints out a message and correlated variable
C ==========================================================
C
C Author: Guenther Dissertori  / 5.1.95
C-----------------------------------------------------------------------

      Implicit NONE
#include "qdecl.h"
#include "qcde.h"

      Real          Var
      Character*(*) Message

      Integer    ErrOut
      Parameter (ErrOut = 33)

#include "qmacro.h"
C-----------------------------------------------------------------------

      write(ErrOut,*)' '
      write(ErrOut,'(X,A,3I8)')'KNEVT,KEVT,KRUN : ',KNEVT,KEVT,KRUN
      write(ErrOut,'(X,A,F10.2)')Message(1:),Var
      write(ErrOut,*)' '

      RETURN
      END



C-----------------------------------------------------------------------
C-----------------------------------------------------------------------



      Subroutine ALEPH_STATUS(IOK,prob_code)
C--------------------------------------------------------------
C
C  Requires  Class 16
C            VDET ITC TPC dE/dx ECAL HCAL to be ON
C 
C  Skip events with readout problems in VDET
C
C  obtained from F.Palla
C--------------------------------------------------------------
      Implicit none
#include "qdecl.h"
#include "qcde.h"

      integer prob_code(10)
      integer i,IOK,jvdet
      real    dumm

#include "qmacro.h"
C--------------------------------------------------------------

      IOK = 1
      do i=1,9
        prob_code(i) = 0 
      enddo
      if(xmcev) then
        IOK = 0
        goto 100
      endif
      if(.not.xvitc)then
        prob_code(1) = 1
        goto 100
      endif
      if(.not.xvtpc)then
        prob_code(2) = 1
        goto 100
      endif
      if(.not.xvtpcd)then
        prob_code(3) = 1
        goto 100
      endif
      if(.not.xvecal)then
        prob_code(4) = 1
        goto 100
      endif
      if(.not.xvhcal)then
        prob_code(5) = 1
        goto 100
      endif
      if(.not.xvdeok(dumm))then
        prob_code(6) = 1
        goto 100
      endif
C
C+ Check Error Status
C      
      if(KEVEES.NE.1)then
        prob_code(7) = 1
        goto 100
      endif
      iok = 0
 100  return
      end



C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      
      SUBROUTINE QCD_Init
C-----------------------------------------------------------------------
C Private initialization step
C ===========================
C
C Track and event selection cuts are initialized. The default values 
C set by this routine can be overwritten with SCUT-cards from the 
C cards file. The final settings of the cuts for the current job are
C reported to the log-file. In addition a set of standard histograms
C for monitoring the the properties of reconstructed tracks and 
C energy-flow objects is booked.
C 
C Author: Michael Schmelling / 21-Jul-1994
C         Guenther Dissertori / 17-Jan-1995 (some histos added)
C
C-----------------------------------------------------------------------

      Implicit NONE
 
      Character*4 Chaint
      Integer     Nlink,Namind
      External    Nlink,Namind,Chaint

      Integer     ind,icut
      Integer     num,kSCUT
      Character*4 cname

C-----------------------------------------------------------------------
C histograms for event monitoring
C ===============================

      Integer idc_nt,idc_d0,idc_z0,idc_th,idc_ph,idc_pt,idc_ee,ide_th
      Integer ide_ph,ide_e3,ide_ee,idh_th,idh_ph,idh_e3,idh_ee
      Integer offs,idev_nch,idev_ech,idev_cts,idev_ntr,offs_ev,offs_evt
      Integer KMCUTS

      Parameter ( idc_nt   =  1 )
      Parameter ( idc_d0   =  2 )
      Parameter ( idc_z0   =  3 )
      Parameter ( idc_th   =  4 )
      Parameter ( idc_ph   =  5 )
      Parameter ( idc_pt   =  6 )
      Parameter ( idc_ee   =  7 )
      Parameter ( ide_th   =  8 )
      Parameter ( ide_ph   =  9 )
      Parameter ( ide_e3   = 10 )
      Parameter ( ide_ee   = 11 )
      Parameter ( idh_th   = 12 )
      Parameter ( idh_ph   = 13 )
      Parameter ( idh_e3   = 14 )
      Parameter ( idh_ee   = 15 )
      Parameter ( offs     = 15 )
      Parameter ( idev_nch = 31 )
      Parameter ( idev_ech = 32 )
      Parameter ( idev_cts = 33 )
      Parameter ( idev_ntr = 34 )
      Parameter ( offs_ev  =  4 )
      Parameter ( offs_evt = 40 )
 
C event class "MC-truth x fiducial cuts"
C ======================================
 
      Parameter ( KMCUTS = -3 )

C Selection cuts
C ==============

      Real            ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
      Common /mscuts/ ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
C-----------------------------------------------------------------------

#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"

C-----------------------------------------------------------------------
C initialize the track and event selection cuts
C =============================================

C--      set default cuts

      CTD0 = 2.     ! maximum d0 for charged tracks 
      CTZ0 = 10.    ! maximum z0 for charged tracks
      CTNT = 4.     ! minimum #tpc hits for a charged track
      CTCT = 0.94   ! maximum abs(cos(th)) of charged tracks 
      CTPT = 0.20   ! minimum pt of charged tracks
      CTCF = 1000.  ! maximum chi2/ndf for track fit
      CCCT = 0.98   ! maximum abs(cos(th)) of ECAL/HCAL objects
      CCEB = 0.40   ! minimum energy ECAL/HCAL objects

      CNCH = 5.     ! minimum number of charged tracks required
      CECH = 15.    ! minimum charged energy required
      CCTS = 0.82   ! maximum abs(cos(th)) of sphericity axis
      CNTR = 13.    ! minimum number of all ch.+neu. tracks

      kSCUT = IW(Namind('SCUT'))
 102  continue
      if(kSCUT.ne.0) then
         cname = Chaint(iw(kSCUT+1))
         if(cname.eq.'CTD0') ctd0 = rw(kSCUT+2)
         if(cname.eq.'CTZ0') ctz0 = rw(kSCUT+2)
         if(cname.eq.'CTNT') ctnt = rw(kSCUT+2)
         if(cname.eq.'CTCT') ctct = rw(kSCUT+2)
         if(cname.eq.'CTPT') ctpt = rw(kSCUT+2)
         if(cname.eq.'CTCF') ctcf = rw(kSCUT+2)
         if(cname.eq.'CCCT') ccct = rw(kSCUT+2)
         if(cname.eq.'CCEB') cceb = rw(kSCUT+2)
         if(cname.eq.'CNCH') cnch = rw(kSCUT+2)
         if(cname.eq.'CECH') cech = rw(kSCUT+2)
         if(cname.eq.'CCTS') ccts = rw(kSCUT+2)
         if(cname.eq.'CNTR') cntr = rw(kSCUT+2)
         kSCUT = iw(kSCUT-1)
         goto 102
      endif 

C--      dump the cut values

      write(6,'(1X,''_QCD_INIT_'',1X,69(''-''))')
      write(6,'( /,1X,A)') 'Selection cuts to be used'
      write(6,'(   1X,A)') '========================='
      write(6,'(   1X,A)') 'criteria for event selection... '
      write(6,'(1X,A,F9.4)') 'Minimum #charged tracks     : ',cnch
      write(6,'(1X,A,F9.4)') 'Minimum charged energy/GeV  : ',cech
      write(6,'(1X,A,F9.4)') 'Maximum abs(cos(th_spher.)) : ',ccts 
      write(6,'(1X,A,F9.4)') 'Minimum #ch+neu  tracks     : ',cntr
      write(6,'(   1X,A)') 'criteria for charged tracks... '
      write(6,'(1X,A,F9.4)') 'Maximum d0                  : ',ctd0
      write(6,'(1X,A,F9.4)') 'Maximum z0                  : ',ctz0
      write(6,'(1X,A,F9.4)') 'Minimum #tpc hits           : ',ctnt
      write(6,'(1X,A,F9.4)') 'Maximum abs(cos(th))        : ',ctct
      write(6,'(1X,A,F9.4)') 'Minimum pt                  : ',ctpt
      write(6,'(1X,A,F9.4)') 'Maxmium Chi2/Ndf            : ',ctcf
      write(6,'(   1X,A)') 'criteria for calorimeter objects... '
      write(6,'(1X,A,F9.4)') 'Maximum abs(cos(th))        : ',ccct
      write(6,'(1X,A,F9.4)') 'Minimum energy in ECAL/HCAL : ',cceb
      write(6,'(1X)')

C-----------------------------------------------------------------------
C book standard monitoring histograms
C ===================================

C     -- all eflow-objects

      call Qbook1(idc_nt,'N(TPC)      CTRK', 22, -0.5, 21.5,0.)
      call Qbook1(idc_d0,'d0          CTRK',100,-10. , 10. ,0.)
      call Qbook1(idc_z0,'z0          CTRK',100,-25. , 25. ,0.)
      call Qbook1(idc_th,'dn/dcos(th) CTRK',180,  0. ,180. ,0.)
      call Qbook1(idc_ph,'dn/dph      CTRK',144,  0. ,360. ,0.)
      call Qbook1(idc_pt,'dn/dpt**2   CTRK',100,  0. ,  1. ,0.)
      call Qbook1(idc_ee,'dn/dE       CTRK',150,  0. , 30. ,0.)

      call Qbook1(ide_th,'dn/dcos(th) ECAL',180,0.,180.,0.)
      call Qbook1(ide_ph,'dn/dph      ECAL',144,0.,360.,0.)
      call Qbook1(ide_e3,'dn/dE**3    ECAL',100,0.,  1.,0.)
      call Qbook1(ide_ee,'dn/dE       ECAL',150,0., 30.,0.)

      call Qbook1(idh_th,'dn/dcos(th) HCAL',180,0.,180.,0.)
      call Qbook1(idh_ph,'dn/dph      HCAL',144,0.,360.,0.)
      call Qbook1(idh_e3,'dn/dE**3    HCAL',100,0.,  1.,0.)
      call Qbook1(idh_ee,'dn/dE       HCAL',100,0., 20.,0.)

C     -- accepted eflow-objects

      call Qbook1(offs+idc_nt,'N(TPC)      accCTRK', 22, -0.5, 21.5,0.)
      call Qbook1(offs+idc_d0,'d0          accCTRK',100,-10. , 10. ,0.)
      call Qbook1(offs+idc_z0,'z0          accCTRK',100,-25. , 25. ,0.)
      call Qbook1(offs+idc_th,'dn/dcos(th) accCTRK',180,  0. ,180. ,0.)
      call Qbook1(offs+idc_ph,'dn/dph      accCTRK',144,  0. ,360. ,0.)
      call Qbook1(offs+idc_pt,'dn/dpt**2   accCTRK',100,  0. ,  1. ,0.)
      call Qbook1(offs+idc_ee,'dn/dE       accCTRK',150,  0. , 30. ,0.)

      call Qbook1(offs+ide_th,'dn/dcos(th) accECAL',180,0.,180.,0.)
      call Qbook1(offs+ide_ph,'dn/dph      accECAL',144,0.,360.,0.)
      call Qbook1(offs+ide_e3,'dn/dE**3    accECAL',100,0.,  1.,0.)
      call Qbook1(offs+ide_ee,'dn/dE       accECAL',150,0., 30.,0.)

      call Qbook1(offs+idh_th,'dn/dcos(th) accHCAL',180,0.,180.,0.)
      call Qbook1(offs+idh_ph,'dn/dph      accHCAL',144,0.,360.,0.)
      call Qbook1(offs+idh_e3,'dn/dE**3    accHCAL',100,0.,  1.,0.)
      call Qbook1(offs+idh_ee,'dn/dE       accHCAL',100,0., 20.,0.)

C     -- accepted eflow-objects for accepted events

      call Qbook1(offs_evt+idc_nt,'N(TPC)      accCTRK acc ev',
     >             22, -0.5, 21.5,0.)
      call Qbook1(offs_evt+idc_d0,'d0          accCTRK acc ev',
     >            100,-10. , 10. ,0.)
      call Qbook1(offs_evt+idc_z0,'z0          accCTRK acc ev',
     >            100,-25. , 25. ,0.)
      call Qbook1(offs_evt+idc_th,'dn/dcos(th) accCTRK acc ev',
     >            180,  0. ,180. ,0.)
      call Qbook1(offs_evt+idc_ph,'dn/dph      accCTRK acc ev',
     >            144,  0. ,360. ,0.)
      call Qbook1(offs_evt+idc_pt,'dn/dpt**2   accCTRK acc ev',
     >            100,  0. ,  1. ,0.)
      call Qbook1(offs_evt+idc_ee,'dn/dE       accCTRK acc ev',
     >            150,  0. , 30. ,0.)

      call Qbook1(offs_evt+ide_th,'dn/dcos(th) accECAL acc ev',
     >            180,0.,180.,0.)
      call Qbook1(offs_evt+ide_ph,'dn/dph      accECAL acc ev',
     >            144,0.,360.,0.)
      call Qbook1(offs_evt+ide_e3,'dn/dE**3    accECAL acc ev',
     >            100,0.,  1.,0.)
      call Qbook1(offs_evt+ide_ee,'dn/dE       accECAL acc ev',
     >            150,0., 30.,0.)

      call Qbook1(offs_evt+idh_th,'dn/dcos(th) accHCAL acc ev',
     >            180,0.,180.,0.)
      call Qbook1(offs_evt+idh_ph,'dn/dph      accHCAL acc ev',
     >            144,0.,360.,0.)
      call Qbook1(offs_evt+idh_e3,'dn/dE**3    accHCAL acc ev',
     >            100,0.,  1.,0.)
      call Qbook1(offs_evt+idh_ee,'dn/dE       accHCAL acc ev',
     >            100,0., 20.,0.)

C     -- histograms for charged-track event properties

      call Qbook1(idev_nch,'EVCUT: nch'        , 60,-0.5, 59.5,0.)
      call Qbook1(idev_ech,'EVCUT: ech'        ,100, 0. ,100. ,0.)
      call Qbook1(idev_cts,'EVCUT: cos(th_sph)', 50, 0. ,  1. ,0.)
      call Qbook1(idev_ntr,'EVCUT: ntrk'       ,100,-0.5, 99.5,0.)

C     -- accepted events

      call Qbook1(offs_ev+idev_nch,'EVCUT: nch, acc ev'
     >            , 60,-0.5,59.5,0.)
      call Qbook1(offs_ev+idev_ech,'EVCUT: ech, acc ev'
     >            ,100,0.,100.,0.)
      call Qbook1(offs_ev+idev_cts,'EVCUT: cos(th_sph), acc ev'
     >            ,50,0.,1.,0.)
      call Qbook1(offs_ev+idev_ntr,'EVCUT: ntrk acc ev'
     >            ,100,-0.5,99.5,0.)

      write(6,'(1X,A,I3,2X,A,/)') 
     >        'Identifiers 1 to',offs_evt+idh_ee,
     >        'reserved for standard track+event selection histograms'
      write(6,'(1X,''_QCD_INIT_'',1X,69(''-''))')
      RETURN
      END 


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE QCD_HadSel(iremc,nch,ech,cts,ntr,passed)
C-----------------------------------------------------------------------
C Check whether the current event is a hadronic event
C ===================================================
C
C Input : I*4 iremc  : KMONTE/KMCUTS/KRECO 
C Output: L*4 passed : .true. if the event passes as a hadronic event
C                    : + info on nch,ech,cts,ntr
C                      = # of charged tracks, energy of charged tracks,
C                        cos(sphericity axis), total number of objects
C
C Description
C -----------
C The information from the E-flow section is used to check whether the 
C current event satisfies the standard charged-track-based hadronic event 
C selection cuts. Depending on the value of IREMC various options are 
C available: IREMC=KMONTE: use all charged stable final state particles
C                  KMCUTS: use all charged tracks within fiducial cuts
C                  KRECO : reconstructed charged tracks within cuts 
C
C Author: Michael Schmelling / 21-Jul-1994
C         Guenther Dissertori / 13-Dec-1994 (some histos added)
C         Guenther Dissertori / 20-Jan-1999 (Output added)
C
C-----------------------------------------------------------------------

      Implicit NONE 

      Integer  iremc
      Logical  passed

      Integer  i
      Integer  nt,idh,nch,nne,ntr
      Real     d0,z0,ct,ph,pt,ee,sph,cts,ech,th,st

C-----------------------------------------------------------------------
C histograms for event monitoring
C ===============================

      Integer idc_nt,idc_d0,idc_z0,idc_th,idc_ph,idc_pt,idc_ee,ide_th
      Integer ide_ph,ide_e3,ide_ee,idh_th,idh_ph,idh_e3,idh_ee
      Integer offs,idev_nch,idev_ech,idev_cts,idev_ntr,offs_ev,offs_evt
      Integer KMCUTS

      Parameter ( idc_nt   =  1 )
      Parameter ( idc_d0   =  2 )
      Parameter ( idc_z0   =  3 )
      Parameter ( idc_th   =  4 )
      Parameter ( idc_ph   =  5 )
      Parameter ( idc_pt   =  6 )
      Parameter ( idc_ee   =  7 )
      Parameter ( ide_th   =  8 )
      Parameter ( ide_ph   =  9 )
      Parameter ( ide_e3   = 10 )
      Parameter ( ide_ee   = 11 )
      Parameter ( idh_th   = 12 )
      Parameter ( idh_ph   = 13 )
      Parameter ( idh_e3   = 14 )
      Parameter ( idh_ee   = 15 )
      Parameter ( offs     = 15 )
      Parameter ( idev_nch = 31 )
      Parameter ( idev_ech = 32 )
      Parameter ( idev_cts = 33 )
      Parameter ( idev_ntr = 34 )
      Parameter ( offs_ev  =  4 )
      Parameter ( offs_evt = 40 )
 
C event class "MC-truth x fiducial cuts"
C ======================================
 
      Parameter ( KMCUTS = -3 )
 
C Selection cuts
C ==============

      Real            ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
      Common /mscuts/ ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
C-----------------------------------------------------------------------

      Logical  xok
      Integer  type

#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"
C-----------------------------------------------------------------------
C branch according to the value of iremc
C ======================================

      passed = .false.


C     -- branch according to the value of iremc

      if(iremc.eq.KMONTE) goto 1111
      if(iremc.eq.KMCUTS) goto 2222
      if(iremc.eq.KRECO ) goto 3333
      
      write(6,'(1X,A,I4,2X,A,I4)') 
     >        'QCD_HADSEL: unknown code IREMC = ',IREMC,
     >                    'treated as KMCUTS = ',KMCUTS  
      goto 2222      

C-----------------------------------------------------------------------
C KMONTE: event selection based on all MC-tracks
C ==============================================

 1111 continue
      nch = 0
      ech = 0.

      call QLZER(KMONTE)

      do i=KFMCT,KLMCT
      if((KSTABC(i).eq. 1) .or.
     >   (KSTABC(i).eq. 2) .or.
     >   (KSTABC(i).eq.-3)     ) then
            if(KCH(i).ne.0) then
               nch = nch + 1
               ech = ech + QE(i)
               call QLTRK(i)
            endif             
         endif
      enddo

C     -- reverse lock and  determine the sphericity axis
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo

      call QJSPHE(sph,'sph_axis',KMONTE)
      if(nch.gt.0) then
         cts = abs(QCT(KPDIR('sph_axis',KMONTE)))    
      else
         cts = 2.
      endif

C -- monitor event properties

      call Hf1(idev_nch,float(nch),1.)
      call Hf1(idev_ech,      ech ,1.)
      call Hf1(idev_cts,      cts ,1.)      

C     -- and decide whether or not the event passes

      if(nch.lt.nint(cnch)) RETURN
      if(ech.lt.     cech ) RETURN
      if(cts.gt.     ccts ) RETURN

      passed = .true.

C -- monitor good events 

      call Hf1(offs_ev+idev_nch,float(nch),1.)
      call Hf1(offs_ev+idev_ech,      ech ,1.)
      call Hf1(offs_ev+idev_cts,      cts ,1.)

      RETURN

C-----------------------------------------------------------------------
C KMCUTS: event selection based on MC-truth x fiducial cuts
C =========================================================

 2222 continue
      nch = 0
      ech = 0.

      call QLZER(KMONTE)

      do i=KFMCT,KLMCT
      if((KSTABC(i).eq. 1) .or.
     >   (KSTABC(i).eq. 2) .or.
     >   (KSTABC(i).eq.-3)     ) then
         if((    KCH(i) .ne.   0) .and.
     >      (abs(QCT(i)).lt.ctct) .and.
     >      (    QPT(i) .gt.ctpt)      ) then
             nch = nch + 1
             ech = ech + QE(i)
             call QLTRK(i)             
         endif
      endif
      enddo

C     -- reverse lock and  determine the sphericity axis
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo  

      call QJSPHE(sph,'sph_axis',KMONTE)
      if(nch.gt.0) then
         cts = abs(QCT(KPDIR('sph_axis',KMONTE)))    
      else
         cts = 2.
      endif

C -- monitor event properties 

      call Hf1(idev_nch,float(nch),1.)
      call Hf1(idev_ech,      ech ,1.)
      call Hf1(idev_cts,      cts ,1.)

C     -- and decide whether or not the event passes

      if(nch.lt.nint(cnch)) RETURN
      if(ech.lt.     cech ) RETURN
      if(cts.gt.     ccts ) RETURN

      passed = .true.

C -- monitor good events 

      call Hf1(offs_ev+idev_nch,float(nch),1.)
      call Hf1(offs_ev+idev_ech,      ech ,1.)
      call Hf1(offs_ev+idev_cts,      cts ,1.)

      RETURN

C-----------------------------------------------------------------------
C KRECO: event selection based on reconstructed information
C =========================================================

 3333 continue

C--       set analysis to charged tracks and open all locks

      call QJOPTR('CH',' ')
      call QLZER(KRECO)

C--       scan over all EFLOW-objects

      nch = 0
      ech = 0.
      nne = 0

      do i=KFEFT,KLEFT
         if(.not.XEFO(i)) then
            write(6,*) 
     >      'QCD_HADSEL: ENFLW-object without EFOL info => ',i
         else
            ct   = QCT(i)
            st   = sqrt(max(1.e-10,1.-ct*ct))
            th   = acos(max(min(ct,1.),-1.)) * 180./QQPI
            ph   = QPH(i) * 180./QQPI
            pt   = QPT(i)
            ee   = QE(i)
            d0   = QDB(i)
            z0   = QZB(i)
            nt   = KFRTNT(i)
            type = KEFOTY(i)

C           -- find good charged tracks, histogram everything
 
            if(     type.lt.0) then
               write(6,*) 
     >         'QCD_EFLOW: unknown EFLOW object type = ',type
            else if(type.le.3) then   ! charged tracks
               call Hf1(idc_d0,      d0 ,1.)
               call Hf1(idc_z0,      z0 ,1.)
               call Hf1(idc_th,      th ,1./st)
               call Hf1(idc_pt,      pt ,1./pt)
               call Hf1(idc_nt,float(nt),1.)
               call Hf1(idc_ph,      ph ,1.)
               call Hf1(idc_ee,      ee ,1.)
               xok = .true.
               if(abs(d0).gt.     ctd0 ) xok = .false.
               if(abs(z0).gt.     ctz0 ) xok = .false.
               if(abs(ct).gt.     ctct ) xok = .false.
               if(    pt .lt.     ctpt ) xok = .false.
               if(    nt .lt.nint(ctnt)) xok = .false.
               if(xok) then 
                  nch = nch + 1
                  ech = ech + ee
                  call QLOCK(i)
               endif
            else if(type.eq.4) then   ! neutral ECAL objects
               xok = .true.
               if(abs(ct).gt.     ccct ) xok = .false.
               if(    ee .lt.     cceb ) xok = .false.
               if(xok) nne = nne + 1
               call Hf1(ide_th,th,1./st)
               call Hf1(ide_ph,ph,1.)
               call Hf1(ide_e3,ee,1./(ee*ee))
               call Hf1(ide_ee,ee,1.)
            else if(type.eq.5) then   ! neutral HCAL objects
               xok = .true.
               if(abs(ct).gt.     ccct ) xok = .false.
               if(    ee .lt.     cceb ) xok = .false.
               if(xok) nne = nne + 1
               call Hf1(idh_th,th,1./st)
               call Hf1(idh_ph,ph,1.)
               call Hf1(idh_e3,ee,1./(ee*ee))
               call Hf1(idh_ee,ee,1.)
            else if(type.gt.8) then
               write(6,*) 
     >         'QCD_HADSEL: unknown EFLOW object type = ',type
            endif
         endif
      enddo

C--        reverse lock and  determine the sphericity axis

  
      call QLREV(KRECO)
      call QJSPHE(sph,'sph_axis',KRECO)
      if(nch.gt.0) then
         cts = abs(QCT(KPDIR('sph_axis',KRECO)))    
      else
         cts = 2.
      endif

      ntr = nch + nne

C--        monitor event properties

      call Hf1(idev_nch,float(nch),1.)
      call Hf1(idev_ech,      ech ,1.)
      call Hf1(idev_cts,      cts ,1.)
      call Hf1(idev_ntr,float(ntr),1.)

C--        decide whether or not it passes...

      if(nch.lt.nint(cnch)) RETURN
      if(ech.lt.     cech ) RETURN
      if(cts.gt.     ccts ) RETURN
      if(ntr.lt.nint(cntr)) RETURN

      passed = .true.

C -- monitor good events 

      call Hf1(offs_ev+idev_nch,float(nch),1.)
      call Hf1(offs_ev+idev_ech,      ech ,1.)
      call Hf1(offs_ev+idev_cts,      cts ,1.)
      call Hf1(offs_ev+idev_ntr,float(ntr),1.)
      RETURN
      END
       
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE QCD_ENFLOW(iremc,csel,ntrk,etot,px,py,pz,passed)
C-----------------------------------------------------------------------
C Select energy-flow objects for further analysis
C ===============================================
C
C Input : I*4   iremc : KRECO/KMCUTS/KMONTE 
C         C*(*) csel  : select 'charged', 'neutral' or 'all' tracks
C         L*4 passed  : .true. if the event passes as a hadronic event
C Output: I*4   ntrk  : number of tracks found for analysis
C         R*4   etot  : total energy 
C         R*4   px    : total x-momentum
C         R*4   py    : total y-momentum
C         R*4   pz    : total z-momentum
C
C The energy flow information is scanned, and all "tracks" which
C pass the selection criteria are locked. Before returning to the 
C calling program the lock-status is reversed. Only the good tracks 
C are unlocked. Routine should be called after QCD_HadSel to be
C sure about meaningful value of passed
C 
C Author: Michael Schmelling / 21-Jul-1994
C         Guenther Dissertori / 17-Jan-1995  ('passed' added)
C-----------------------------------------------------------------------

      Implicit NONE

      Integer       iremc,ntrk
      Character*(*) csel
      Real          etot,px,py,pz
 
      Integer       i

      Logical       xok,passed
      Integer       nt,type,isel
      Real          ct,ph,d0,z0,ee,pt,st,th

C-----------------------------------------------------------------------
C histograms for event monitoring
C ===============================

      Integer idc_nt,idc_d0,idc_z0,idc_th,idc_ph,idc_pt,idc_ee,ide_th
      Integer ide_ph,ide_e3,ide_ee,idh_th,idh_ph,idh_e3,idh_ee
      Integer offs,idev_nch,idev_ech,idev_cts,idev_ntr,offs_ev,offs_evt
      Integer KMCUTS

      Parameter ( idc_nt   =  1 )
      Parameter ( idc_d0   =  2 )
      Parameter ( idc_z0   =  3 )
      Parameter ( idc_th   =  4 )
      Parameter ( idc_ph   =  5 )
      Parameter ( idc_pt   =  6 )
      Parameter ( idc_ee   =  7 )
      Parameter ( ide_th   =  8 )
      Parameter ( ide_ph   =  9 )
      Parameter ( ide_e3   = 10 )
      Parameter ( ide_ee   = 11 )
      Parameter ( idh_th   = 12 )
      Parameter ( idh_ph   = 13 )
      Parameter ( idh_e3   = 14 )
      Parameter ( idh_ee   = 15 )
      Parameter ( offs     = 15 )
      Parameter ( idev_nch = 31 )
      Parameter ( idev_ech = 32 )
      Parameter ( idev_cts = 33 )
      Parameter ( idev_ntr = 34 )
      Parameter ( offs_ev  =  4 )
      Parameter ( offs_evt = 40 )
 
C event class "MC-truth x fiducial cuts"
C ======================================
 
      Parameter ( KMCUTS = -3 )

C Selection cuts
C ==============

      Real            ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
      Common /mscuts/ ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
C-----------------------------------------------------------------------

#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"

C-----------------------------------------------------------------------
C Initialization
C ==============
 
C     -- global properties from all accepted tracks

      ntrk = 0
      etot = 0.
      px   = 0.
      py   = 0.
      pz   = 0.

C     -- decode csel

      isel = 0
      if(csel(1:1).eq.'c'.or.csel(1:1).eq.'C') isel = 1
      if(csel(1:1).eq.'n'.or.csel(1:1).eq.'N') isel = 2
      if(csel(1:1).eq.'a'.or.csel(1:1).eq.'A') isel = 3
      if(isel.eq.0) then
         isel = 3
         write(6,'(1X,A,A,2X,A)') 
     >           'QCD_ENFLOW: unknown track selection ',csel,
     >           'assumed to be ALL' 
      endif

C     -- branch according to the value of iremc

      if(iremc.eq.KMONTE) goto 1111
      if(iremc.eq.KMCUTS) goto 2222
      if(iremc.eq.KRECO ) goto 3333
      
      write(6,'(1X,A,I4,2X,A,I4)') 
     >        'QCD_EFLOW: unknown code IREMC = ',IREMC,
     >                    'treated as KMCUTS = ',KMCUTS  
      goto 2222      

C-----------------------------------------------------------------------
C pure MC-truth (based on KSTABC=1,2,-3)
C =============

 1111 continue

      call QLZER(KMONTE)

C     -- loop over all MC-truth tracks 

      do i=KFMCT,KLMCT
         xok = .false.
         if((KSTABC(i).eq. 1) .or.
     >      (KSTABC(i).eq. 2) .or.
     >      (KSTABC(i).eq.-3)     ) then
            xok = .true.
            if(KCH(i).ne.0.and.isel.eq.2) xok = .false.
            if(KCH(i).eq.0.and.isel.eq.1) xok = .false.      
         endif

C        -- sum up global properties of accepted Eflow-objects

         if(xok) then
            call QLTRK(i)
            ntrk = ntrk + 1
            etot = etot + QE(i)
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
         endif

      enddo

C--      reverse lock status and return
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo

      RETURN

C-----------------------------------------------------------------------
C fiducial cuts EFLOW for MC-truth (based on KSTABC=1,2,-3)
C ================================

 2222 continue

      call QLZER(KMONTE)

C     -- loop over all MC-truth tracks 

      do i=KFMCT,KLMCT
         xok = .false.
         if((KSTABC(i).eq. 1) .or.
     >      (KSTABC(i).eq. 2) .or.
     >      (KSTABC(i).eq.-3)     ) then

            ct = abs(QCT(i))
            pt = QPT(i)
            ee = QE(i)

C           -- check: invisible,sical,lcal,charged,cal-object

            if((Index(CQTPN(i),'nu').eq.1).or.(ct.gt.ccct)) then
               xok = .false.
            else 
               if(KCH(i).ne.0) then 
                  if((pt.gt.ctpt).and.(ct.lt.ctct)) xok = .true.
               else
                  if((ee.gt.cceb).and.(ct.lt.ccct)) xok = .true.
               endif
            endif
C           -----
            if(KCH(i).ne.0.and.isel.eq.2) xok = .false.
            if(KCH(i).eq.0.and.isel.eq.1) xok = .false.      
         endif

C        -- sum up global properties of accepted Eflow-objects

         if(xok) then
            call QLTRK(i)
            ntrk = ntrk + 1
            etot = etot + ee
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
         endif

      enddo

C--      reverse lock status and return
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo

      RETURN

C-----------------------------------------------------------------------
C scan EFLOW-objects from reconstructed information
C =================================================

 3333 continue

C     -- reset all locks

      call QLZER(KRECO)

      do i=KFEFT,KLEFT
         xok = .false.
         if(.not.XEFO(i)) then
            write(6,*) 'QCD_EFLOW: ENFLW-object without EFOL info => ',i
         else
            ct   = QCT(i)
            st   = sqrt(max(1.e-10,1.-ct*ct))
            th   = acos(max(min(ct,1.),-1.)) * 180./QQPI
            ph   = QPH(i) * 180./QQPI
            pt   = QPT(i)
            ee   = QE(i)
            d0   = QDB(i)
            z0   = QZB(i)
            nt   = KFRTNT(i)
            type = KEFOTY(i)

C           -- check cuts
 
            if(     type.lt.0) then
               write(6,*) 'QCD_EFLOW: unknown EFLOW object type = ',type
            else if(type.le.3) then
               xok = .true.
               if(abs(d0).gt.     ctd0 ) xok = .false.
               if(abs(z0).gt.     ctz0 ) xok = .false.
               if(abs(ct).gt.     ctct ) xok = .false.
               if(abs(pt).lt.     ctpt ) xok = .false.
               if(    nt .lt.nint(ctnt)) xok = .false.
            else if(type.le.5) then
               xok = .true.
               if(abs(ct).gt.     ccct ) xok = .false.
               if(    ee .lt.     cceb ) xok = .false.
            else if(type.gt.8) then
               write(6,*) 'QCD_EFLOW: unknown EFLOW object type = ',type
            endif
         endif

         if(type.le.3.and.isel.eq.2) xok = .false.
         if(type.eq.4.and.isel.eq.1) xok = .false.
         if(type.eq.5.and.isel.eq.1) xok = .false.

C        -- accumulate global event properties from good eflow-objects

         if(xok) then
            call QLOCK(i)
            ntrk = ntrk + 1
            etot = etot + ee
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
            if(     type.le.3) then
               call Hf1(offs+idc_d0,d0,1.)
               call Hf1(offs+idc_z0,z0,1.)
               call Hf1(offs+idc_th,th,1./st)
               call Hf1(offs+idc_pt,pt,1./pt)
               call Hf1(offs+idc_nt,float(nt),1.)
               call Hf1(offs+idc_ph,ph,1.)
               call Hf1(offs+idc_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+idc_d0,d0,1.)
                 call Hf1(offs_evt+idc_z0,z0,1.)
                 call Hf1(offs_evt+idc_th,th,1./st)
                 call Hf1(offs_evt+idc_pt,pt,1./pt)
                 call Hf1(offs_evt+idc_nt,float(nt),1.)
                 call Hf1(offs_evt+idc_ph,ph,1.)
                 call Hf1(offs_evt+idc_ee,ee,1.)
               endif
            else if(type.eq.4) then
               call Hf1(offs+ide_th,th,1./st)
               call Hf1(offs+ide_ph,ph,1.)
               call Hf1(offs+ide_e3,ee,1./(ee*ee))
               call Hf1(offs+ide_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+ide_th,th,1./st)
                 call Hf1(offs_evt+ide_ph,ph,1.)
                 call Hf1(offs_evt+ide_e3,ee,1./(ee*ee))
                 call Hf1(offs_evt+ide_ee,ee,1.)
               endif
            else if(type.eq.5) then
               call Hf1(offs+idh_th,th,1./st)
               call Hf1(offs+idh_ph,ph,1.)
               call Hf1(offs+idh_e3,ee,1./(ee*ee))
               call Hf1(offs+idh_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+idh_th,th,1./st)
                 call Hf1(offs_evt+idh_ph,ph,1.)
                 call Hf1(offs_evt+idh_e3,ee,1./(ee*ee))
                 call Hf1(offs_evt+idh_ee,ee,1.)                 
               endif
            endif
         endif

      enddo

C--      reverse lock status and return
      
      call QLREV(KRECO)
      RETURN

      END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE GD_ENFLOW(iremc,csel,ntrk,etot,px,py,pz,passed)
C-----------------------------------------------------------------------
C Select energy-flow objects for further analysis
C ===============================================
C
C  Almost the same as QCD_ENFLW, but in case of KMONTE the
C  parton level can be selected!!
C
C Input : I*4   iremc : KRECO/KMCUTS/KMONTE 
C         C*(*) csel  : select 'charged', 'neutral' , 'all' tracks or 'parton' level
C         L*4 passed  : .true. if the event passes as a hadronic event
C Output: I*4   ntrk  : number of tracks found for analysis
C         R*4   etot  : total energy 
C         R*4   px    : total x-momentum
C         R*4   py    : total y-momentum
C         R*4   pz    : total z-momentum
C
C The energy flow information is scanned, and all "tracks" which
C pass the selection criteria are locked. Before returning to the 
C calling program the lock-status is reversed. Only the good tracks 
C are unlocked. Routine should be called after QCD_HadSel to be
C sure about meaningful value of passed
C 
C Author: Michael Schmelling / 21-Jul-1994
C         Guenther Dissertori / 17-Jan-1995  ('passed' added)
C-----------------------------------------------------------------------

      Implicit NONE

      Integer       iremc,ntrk
      Character*(*) csel
      Real          etot,px,py,pz
 
      Integer       i,Istring

      Logical       xok,passed,IsParton
      Integer       nt,type,isel
      Real          ct,ph,d0,z0,ee,pt,st,th

C-----------------------------------------------------------------------
C histograms for event monitoring
C ===============================

      Integer idc_nt,idc_d0,idc_z0,idc_th,idc_ph,idc_pt,idc_ee,ide_th
      Integer ide_ph,ide_e3,ide_ee,idh_th,idh_ph,idh_e3,idh_ee
      Integer offs,idev_nch,idev_ech,idev_cts,idev_ntr,offs_ev,offs_evt
      Integer KMCUTS

      Parameter ( idc_nt   =  1 )
      Parameter ( idc_d0   =  2 )
      Parameter ( idc_z0   =  3 )
      Parameter ( idc_th   =  4 )
      Parameter ( idc_ph   =  5 )
      Parameter ( idc_pt   =  6 )
      Parameter ( idc_ee   =  7 )
      Parameter ( ide_th   =  8 )
      Parameter ( ide_ph   =  9 )
      Parameter ( ide_e3   = 10 )
      Parameter ( ide_ee   = 11 )
      Parameter ( idh_th   = 12 )
      Parameter ( idh_ph   = 13 )
      Parameter ( idh_e3   = 14 )
      Parameter ( idh_ee   = 15 )
      Parameter ( offs     = 15 )
      Parameter ( idev_nch = 31 )
      Parameter ( idev_ech = 32 )
      Parameter ( idev_cts = 33 )
      Parameter ( idev_ntr = 34 )
      Parameter ( offs_ev  =  4 )
      Parameter ( offs_evt = 40 )
 
C event class "MC-truth x fiducial cuts"
C ======================================
 
      Parameter ( KMCUTS = -3 )

C Selection cuts
C ==============

      Real            ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
      Common /mscuts/ ctd0,ctz0,ctnt,ctct,ctpt,ctcf,
     >                ccct,cceb,
     >                cnch,cech,ccts,cntr
C-----------------------------------------------------------------------

#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"

C-----------------------------------------------------------------------
C Initialization
C ==============
 
C     -- global properties from all accepted tracks

      ntrk = 0
      etot = 0.
      px   = 0.
      py   = 0.
      pz   = 0.

C     -- decode csel

      isel = 0
      if(csel(1:1).eq.'c'.or.csel(1:1).eq.'C') isel = 1
      if(csel(1:1).eq.'n'.or.csel(1:1).eq.'N') isel = 2
      if(csel(1:1).eq.'a'.or.csel(1:1).eq.'A') isel = 3
      if(csel(1:1).eq.'p'.or.csel(1:1).eq.'P') isel = 4
      if(isel.eq.0) then
         isel = 3
         write(6,'(1X,A,A,2X,A)') 
     >           'GD_ENFLOW: unknown track selection ',csel,
     >           'assumed to be ALL' 
      endif

C     -- branch according to the value of iremc

      if(iremc.eq.KMONTE) goto 1111
      if(iremc.eq.KMCUTS) goto 2222
      if(iremc.eq.KRECO ) goto 3333
      
      write(6,'(1X,A,I4,2X,A,I4)') 
     >        'GD_EFLOW: unknown code IREMC = ',IREMC,
     >                    'treated as KMCUTS = ',KMCUTS  
      goto 2222      

C-----------------------------------------------------------------------
C pure MC-truth (based on KSTABC=1,2,-3)
C =============

 1111 continue

      call QLZER(KMONTE)

C     -- loop over all MC-truth tracks 

C -- if parton level is requested
      if(isel.eq.4) then

C -- find string 
        do i=KFMCT,KLMCT
          if (KTPCOD(i).eq.380.and.CQTPN(i).eq.'string') then
            Istring = i
            goto 10
          endif
        enddo
 10     continue
C -- find parton level
        do i=KFMCT,KLMCT       
          IsParton = .False.
          IsParton = ((i.lt.Istring).and.
     &                (KLUNDS(i).eq.11.or.KLUNDS(i).eq.12).and.
     &                (KTPCOD(i).ge.100.and.KTPCOD(i).lt.112)).or.
C -- check especially for photons
     &               ((i.lt.Istring).and.
     &                (KLUNDS(i).eq.1).and.(KTPCOD(i).eq.1))

C        -- sum up global properties of accepted Eflow-objects

          if(IsParton) then
            call QLTRK(i)
            ntrk = ntrk + 1
            etot = etot + QE(i)
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
          endif
        enddo

      else			
C -- else select stable stuff
        do i=KFMCT,KLMCT
          xok = .false.
          if((KSTABC(i).eq. 1) .or.
     >         (KSTABC(i).eq. 2) .or.
     >         (KSTABC(i).eq.-3)     ) then
            xok = .true.
            if(KCH(i).ne.0.and.isel.eq.2) xok = .false.
            if(KCH(i).eq.0.and.isel.eq.1) xok = .false.      
          endif

C        -- sum up global properties of accepted Eflow-objects

          if(xok) then
            call QLTRK(i)
            ntrk = ntrk + 1
            etot = etot + QE(i)
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
          endif
        enddo
      endif



C--      reverse lock status and return
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo

      RETURN

C-----------------------------------------------------------------------
C fiducial cuts EFLOW for MC-truth (based on KSTABC=1,2,-3)
C ================================

 2222 continue

      call QLZER(KMONTE)

C     -- loop over all MC-truth tracks 

      do i=KFMCT,KLMCT
         xok = .false.
         if((KSTABC(i).eq. 1) .or.
     >      (KSTABC(i).eq. 2) .or.
     >      (KSTABC(i).eq.-3)     ) then

            ct = abs(QCT(i))
            pt = QPT(i)
            ee = QE(i)

C           -- check: invisible,sical,lcal,charged,cal-object

            if((Index(CQTPN(i),'nu').eq.1).or.(ct.gt.ccct)) then
               xok = .false.
            else 
               if(KCH(i).ne.0) then 
                  if((pt.gt.ctpt).and.(ct.lt.ctct)) xok = .true.
               else
                  if((ee.gt.cceb).and.(ct.lt.ccct)) xok = .true.
               endif
            endif
C           -----
            if(KCH(i).ne.0.and.isel.eq.2) xok = .false.
            if(KCH(i).eq.0.and.isel.eq.1) xok = .false.      
         endif

C        -- sum up global properties of accepted Eflow-objects

         if(xok) then
            call QLTRK(i)
            ntrk = ntrk + 1
            etot = etot + ee
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
         endif

      enddo

C--      reverse lock status and return
c      call QLREV(KMONTE)
      do i=KFMCT,KLMCT
        if (XLOCK(i)) then
          call QLUTRK(i)
        else
          call QLTRK(i)
        endif
      enddo

      RETURN

C-----------------------------------------------------------------------
C scan EFLOW-objects from reconstructed information
C =================================================

 3333 continue

C     -- reset all locks

      call QLZER(KRECO)

      do i=KFEFT,KLEFT
         xok = .false.
         if(.not.XEFO(i)) then
            write(6,*) 'GD_EFLOW: ENFLW-object without EFOL info => ',i
         else
            ct   = QCT(i)
            st   = sqrt(max(1.e-10,1.-ct*ct))
            th   = acos(max(min(ct,1.),-1.)) * 180./QQPI
            ph   = QPH(i) * 180./QQPI
            pt   = QPT(i)
            ee   = QE(i)
            d0   = QDB(i)
            z0   = QZB(i)
            nt   = KFRTNT(i)
            type = KEFOTY(i)

C           -- check cuts
 
            if(     type.lt.0) then
               write(6,*) 'GD_EFLOW: unknown EFLOW object type = ',type
            else if(type.le.3) then
               xok = .true.
               if(abs(d0).gt.     ctd0 ) xok = .false.
               if(abs(z0).gt.     ctz0 ) xok = .false.
               if(abs(ct).gt.     ctct ) xok = .false.
               if(abs(pt).lt.     ctpt ) xok = .false.
               if(    nt .lt.nint(ctnt)) xok = .false.
            else if(type.le.5) then
               xok = .true.
               if(abs(ct).gt.     ccct ) xok = .false.
               if(    ee .lt.     cceb ) xok = .false.
            else if(type.gt.8) then
               write(6,*) 'GD_EFLOW: unknown EFLOW object type = ',type
            endif
         endif

         if(type.le.3.and.isel.eq.2) xok = .false.
         if(type.eq.4.and.isel.eq.1) xok = .false.
         if(type.eq.5.and.isel.eq.1) xok = .false.

C        -- accumulate global event properties from good eflow-objects

         if(xok) then
            call QLOCK(i)
            ntrk = ntrk + 1
            etot = etot + ee
            px   = px   + QX(i)
            py   = py   + QY(i)
            pz   = pz   + QZ(i)
            if(     type.le.3) then
               call Hf1(offs+idc_d0,d0,1.)
               call Hf1(offs+idc_z0,z0,1.)
               call Hf1(offs+idc_th,th,1./st)
               call Hf1(offs+idc_pt,pt,1./pt)
               call Hf1(offs+idc_nt,float(nt),1.)
               call Hf1(offs+idc_ph,ph,1.)
               call Hf1(offs+idc_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+idc_d0,d0,1.)
                 call Hf1(offs_evt+idc_z0,z0,1.)
                 call Hf1(offs_evt+idc_th,th,1./st)
                 call Hf1(offs_evt+idc_pt,pt,1./pt)
                 call Hf1(offs_evt+idc_nt,float(nt),1.)
                 call Hf1(offs_evt+idc_ph,ph,1.)
                 call Hf1(offs_evt+idc_ee,ee,1.)
               endif
            else if(type.eq.4) then
               call Hf1(offs+ide_th,th,1./st)
               call Hf1(offs+ide_ph,ph,1.)
               call Hf1(offs+ide_e3,ee,1./(ee*ee))
               call Hf1(offs+ide_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+ide_th,th,1./st)
                 call Hf1(offs_evt+ide_ph,ph,1.)
                 call Hf1(offs_evt+ide_e3,ee,1./(ee*ee))
                 call Hf1(offs_evt+ide_ee,ee,1.)
               endif
            else if(type.eq.5) then
               call Hf1(offs+idh_th,th,1./st)
               call Hf1(offs+idh_ph,ph,1.)
               call Hf1(offs+idh_e3,ee,1./(ee*ee))
               call Hf1(offs+idh_ee,ee,1.)
               if (passed) then
                 call Hf1(offs_evt+idh_th,th,1./st)
                 call Hf1(offs_evt+idh_ph,ph,1.)
                 call Hf1(offs_evt+idh_e3,ee,1./(ee*ee))
                 call Hf1(offs_evt+idh_ee,ee,1.)                 
               endif
            endif
         endif

      enddo

C--      reverse lock status and return
      
      call QLREV(KRECO)
      RETURN

      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE QCD_EvDump
C-----------------------------------------------------------------------
C Dump the current event
C ======================
C
C Author: Michael Schmelling / 17-Dec-1993
C-----------------------------------------------------------------------

      Implicit NONE

      Integer  i,imoth

#include "qdecl.h"
#include "qcde.h"
#include "qmacro.h"
C-----------------------------------------------------------------------

      write(6,'(//,1X,A,2I10)') 'Event dump: run, event: ',KRUN,KEVT

C--      dump MC-truth information

      if(XMCEV) write(6,'(/,1X,A,/)') 'MC truth information'
      do i=KFMCT,KLMCT
         imoth = 0
         if(KNMOTH(i).gt.0) imoth = KMOTH(i,1)
         write(6,'(1X,A,2I6,2X,A12,2X,A,4F10.4,2X,I6)')
     >           'track,mother: ',i,imoth,CQTPN(i), 
     >           'px,py,pz,e,klund: ',QX(i),QY(i),QZ(i),QE(i),KLUNDS(i)
      enddo

C--      dump charged track information

      write(6,'(//,1X,A,/)') 'Reconstructed charged tracks'
      do i=KFCHT,KLCHT
         write(6,'(1X,A,I6,2X,A,2F10.4,2X,A,4F10.4,I6)')
     >           'track: ',i,
     >           'd0,z0: ',QDB(i),QZB(i),
     >           'px,py,pz,e,q: ',QX(i),QY(i),QZ(i),QE(i),KCH(i)
      enddo    
  
C--      dump eflow information

      write(6,'(//,1X,A,/)') 'Reconstructed energy flow objects'
      do i=KFEFT,KLEFT
         write(6,'(1X,A,I6,2X,A,2F10.4,2X,A,4F10.4,I6,L4)')
     >           'track: ',i,
     >           'd0,z0: ',QDB(i),QZB(i),
     >           'px,py,pz,e,type,lock: ',
     >            QX(i),QY(i),QZ(i),QE(i),KEFOTY(i),XLOCK(i)
      enddo    

      write(6,'(1X,79(''-''))')

      RETURN
      END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE Calc_Evshapes(iret)
C-----------------------------------------------------------------------
C ALPHA user analysis of event shapes and jet rates
C =================================================
C
C IO : iret on output : 0 = OK
C
C Author: Guenther Dissertori  / 28.5.2001 
C-----------------------------------------------------------------------

      Implicit NONE

#include "qdecl.h"
#include "qcde.h"

      Integer   iret

      Integer   NefoMax
      Parameter (NefoMax = 300)

      Logical   Miss
      Integer   ntrk,itkl,i,j,IretGlobal
      Integer   Scheme,Error

      Real      Pefo(4,NefoMax),Pjet(4,NefoMax),Ecm
      Real      EvS(20),ThD(4,3),Yn(NefoMax),xth
      Real      Abz,Anr,Aks,A34
      
      Integer   Nhad
      Common /Norm/ Nhad

      Integer   ncut,icut,Njet,iddp
      Parameter(ncut=3)      
      Logical  Is3jet(ncut)
      Real     yc

C -- histo definitions
      Include 'histonaming.h'

#include "qmacro.h"
C-----------------------------------------------------------------------

C -- initialization
      iret = 1
      Ecm = QELEP

C -- Copy momenta to internal array
      itkl = 0
      do i=KFEFT,KLEFT
        if(.not.XEFO(i)) then
          call ErrMess('Calc_Evshapes: ENFLW-object without EFOL info,
     >Index of object i=',float(i))
        else if (.not.XLOCK(i)) then
          itkl = itkl+1

          Pefo(1,itkl) = QX(i)
          Pefo(2,itkl) = QY(i)
          Pefo(3,itkl) = QZ(i)
          Pefo(4,itkl) = QE(i)

C -- protection against abnormally large events
          if (itkl.eq.NefoMax) then
            call ErrMess('Calc_Evshapes, limit NefoMax reached!',
     >           float(itkl))
            IretGlobal = -15
            call HFILL(idc1,float(IretGlobal),0.,1.)
            RETURN  
          endif

        endif
      enddo
      ntrk = itkl
 
C-- do jet clustering
      Scheme = 1         ! Durham E-scheme
      do icut=1,ncut
        Is3jet(icut) = .False.
        yc = Ycut3j(icut)
        call DurhamJn(Pefo,ntrk,Scheme,yc,Njet,Pjet,Error)
        if (Error.ne.0) then
          call ErrMess('Calc_Evshapes, Error in clustering=',
     >         float(Error))
          IretGlobal = -7
          call HFILL(idc1,float(IretGlobal),0.,1.)
          RETURN  
        else
          if (Njet.eq.3) then
            Is3jet(icut) = .True.
            call HFILL(idjr,float(icut),0.,1.)
          endif
        endif
      enddo


C -- do not consider missing momentum
      Miss = .False.

C -- Calculate Eventshapes for all particles
      call EvShapes(Pefo,ntrk,Miss,EvS,ThD,Error)

      if (Error.ne.0) then
        call ErrMess('Calc_EvShapes: error occured,Error=',
     >                float(Error))
        IretGlobal = -8
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN   ! ------>
      endif

C -- for -ln(y3) we use this routine
      call DurhamYn(Pefo,ntrk,Scheme,Yn,Error)
      if (Error.ne.0) then
        call ErrMess('Calc_Evshapes, Error in Yn clustering=',
     >       float(Error))
        IretGlobal = -9
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                  ! ----->
      endif
      if (Yn(1).le.0.) then
        call ErrMess('CalcEvshapes, strange Y3 value=',Yn(1))
        IretGlobal = -10
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                  ! ----->
      endif

C -- bookkeeping of events for which eventshape calc. is ok
      call HFILL(idc2,3.,0.,1.)

C -- fill the histos 
C -- see routine EvShapes for definition of the array EvS
      call HFILL(idth,EvS(1),0.,1.)
      call HFILL(idta,EvS(2),0.,1.)
      call HFILL(idtm,EvS(3),0.,1.)
      call HFILL(idmh,EvS(6),0.,1.)
      call HFILL(idcp,EvS(5),0.,1.)
      call HFILL(idbt,EvS(9),0.,1.)
      call HFILL(idbw,EvS(10),0.,1.)
      call HFILL(idly,-log(Yn(1)),0.,1.)

C -- fill the Dparameter values
      call HFILL(iddi,EvS(11),0.,1.)
      do icut=1,ncut
        if     (icut.eq.1) then
          iddp = idd1
        elseif (icut.eq.2) then
          iddp = idd2
        elseif (icut.eq.3) then
          iddp = idd3
        endif
        if (Is3jet(icut)) then 
          call HFILL(iddp,EvS(11),0.,1.)
        endif
      enddo

C -- fill the moments of 1-Thrust
      xth = 1-EvS(1)
      do i=1,4
        call HFILL(idmo,float(i),0.,xth**i)
        call HFILL(idmo,float(i+4),0.,xth**(2*i))
      enddo

C -- angular correlations in 4-jet events
C -- check first if it is a 4-jet event
      Scheme = 1         ! Durham E-scheme
      yc = Ycut4j
      call DurhamJn(Pefo,ntrk,Scheme,yc,Njet,Pjet,Error)
      if (Error.ne.0) then
        call ErrMess('Calc_Evshapes, Error in 4-jet clustering=',
     >       float(Error))
        IretGlobal = -11
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN   ! ----->
      else
C -- if it is a 4-jet event, compute the angles and fill them
        if (Njet.eq.4) then
          call Fourjets(Pjet,Abz,Aks,Anr,A34,Error)
          if (Error.ne.0) then
            call ErrMess('Calc_Evshapes, Error in 4-jet angles=',
     >       float(Error))
            IretGlobal = -12
            call HFILL(idc1,float(IretGlobal),0.,1.)
            RETURN              ! ----->
          else
            call HFILL(idbz,abs(cos(Abz)),0.,1.)
            call HFILL(idnr,abs(cos(Anr)),0.,1.)
            call HFILL(idks,    cos(Aks) ,0.,1.)
            call HFILL(id34,    cos(A34) ,0.,1.)

C -- bookkeeping of events for which 4-jet angular corr. are OK
            call HFILL(idc2,4.,0.,1.)

          endif
        endif
      endif

C -- OK, everthing done
      iret = 0

 9999 RETURN 
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


      SUBROUTINE Calc_Inclusives(iret)
C-----------------------------------------------------------------------
C ALPHA user analysis of inclusive charged particle distributions
C ===============================================================
C
C IO : iret on output : 0 = OK
C
C Author: Guenther Dissertori  / 28.5.2001 
C-----------------------------------------------------------------------

      Implicit NONE

#include "qdecl.h"
#include "qcde.h"

      Integer   iret,IretGlobal

      Real      Ecm,eigval,THRUST
      Real      xp,xsi,ptth,yra,pti,pto,xup,xdo
      Integer   i,imajo,imidl,imino,ITHRU

      Integer   Nhad
      Common /Norm/ Nhad

C -- histo definitions
      Include 'histonaming.h'

#include "qmacro.h"
C-----------------------------------------------------------------------

C -- initialization
      iret = 1
      Ecm = QELEP

C -- compute event plane axes (from Sphericity tensor) and Thrust, 
C -- with charged tracks only !!
      call QJOPTR('CH',' ')
      call QJEIG (eigval,'eigvec',KRECO)
      imajo = KPDIR ('eigvec',KRECO)
      imidl = KFOLLO (imajo)
      imino = KFOLLO (imidl)

      call QJTHRU(THRUST,'THRUSTCHG',KRECO)
      if (THRUST.eq.0.) then
        IretGlobal = -13
        call HFILL(idc1,float(IretGlobal),0.,1.)
        write(6,*) 'Calc_Inclusives : No thrust',KRUN,KEVT
        RETURN                ! ------------>
      endif
      ITHRU = KPDIR('THRUSTCHG',KRECO)
      if (ITHRU.le.0.or.QP(ITHRU).eq.0.) THEN
        CALL QWMESE('Calc_Inclusives:Problems with thrust axis')
        IretGlobal = -14
        call HFILL(idc1,float(IretGlobal),0.,1.)
        RETURN                ! ------------>
      endif
C -- store the angle of the axis with charged tracks only
      call HFILL(idt3,abs(QCT(ITHRU)),0.,1.)

C -- now the loop over charged eflow objects
      do i=KFEFT,KLEFT
        if(.not.XEFO(i)) then
          call ErrMess('Calc_Inclusives: ENFLW-object without EFOL info,
     >Index of object i=',float(i))
C -- check if a good one
        else if (.not.XLOCK(i)) then

C -- if charged
          if (KEFOTY(i).le.3) then

            xp    = 2*QP(i)/Ecm
            xsi   = -log(xp)
            
C -- momentum component parallel to thrust
            ptth  = abs(QPPAR(i,ITHRU))
C -- rapidity w.r.t. thrust axis
            xup   = QE(i) + ptth
            xdo   = QE(i) - ptth
            if (xdo.eq.0..or.xup.eq.0.) then
              yra = -1.
            else
              yra   = 0.5*log(xup/xdo)
            endif
C -- momenta w.r.t. event plane
            pti = abs(QPPAR(i,imajo))
            pto = abs(QPPAR(i,imino))

C -- fill the histos
            call HFILL(idxp,xp ,0.,1.)
            call HFILL(idxi,xsi,0.,1.)
            call HFILL(idra,yra,0.,1.)
            call HFILL(idpi,pti,0.,1.)
            call HFILL(idpo,pto,0.,1.)
          endif
        endif
      enddo


C -- OK, everthing done
      iret = 0

 9999 RETURN 
      END
      

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


C include additional code : the utilities package
C =======================
      Include 'qcd_utilities.F'


