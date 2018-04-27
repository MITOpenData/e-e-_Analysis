C-------------------------------------------------
C include file which defines naming and numbering
C of historgrams used in the QCD example program
C
C Author : G. Dissertori,  23-MAY-01
C-------------------------------------------------

C -- event shapes
      Integer idth,
     >        idta,
     >        idtm,
     >        idmh,
     >        idcp,
     >        idbt,
     >        idbw,
     >        idly
      Integer offevs

      Parameter (idth = 100)          ! Thrust
      Parameter (idta = 110)          ! Thrust Major   
      Parameter (idtm = 120)          ! Thrust Minor   
      Parameter (idmh = 130)          ! Heavy jet mass
      Parameter (idcp = 140)          ! Cparameter   
      Parameter (idbt = 150)          ! Total jet broadening
      Parameter (idbw = 160)          ! Wide jet broadening
      Parameter (idly = 170)          ! -ln(y3) : differential 2-jet rate
      Parameter (offevs = 0)          ! offset

C -- moments of event shapes
      Integer idmo
      Integer offmom

      Parameter (idmo = 200)          ! moments of Thrust
      Parameter (offmom = 0)          ! offset

C -- jet rates
      Integer idjr
      Integer offjra
      Real    Ycut3j(3)

      Parameter (idjr = 300)          ! 3-jet rate for three ycut values
      Parameter (offjra = 0)          ! offset
      Data Ycut3j / 0.1, 0.05, 0.025 /! the ycut values

C -- D parameter for different topologies 
      Integer iddi,
     >        idd1,
     >        idd2,
     >        idd3
      Integer offdpa

      Parameter (iddi = 180)          ! Dparameter, inclusively
      Parameter (idd1 = 181)          ! for ycut = Ycut3j(1)
      Parameter (idd2 = 182)          !                   2
      Parameter (idd3 = 183)          !                   3
      Parameter (offdpa = 0)          ! offset

C -- Angular correlations in 4-jet events
      Integer idbz,
     >        idnr,
     >        idks,
     >        id34
      Integer offang
      Real Ycut4j                     ! ycut for 4-jet events

      Parameter (idbz = 410)          ! Bengtsson-Zerwas angle
      Parameter (idnr = 420)          ! Nachtmann-Reiter
      Parameter (idks = 430)          ! Koerner-Schierholz-Willrodt
      Parameter (id34 = 440)          ! alpha_34
      Parameter (offang = 0)          ! offset
      Data Ycut4j / 0.008 /           ! the ycut value

C -- inclusive charged particle distributions
      Integer idxp,
     >        idxi,
     >        idra,
     >        idpi,
     >        idpo
      Integer offinc

      Parameter (idxp = 510)          ! xp=pch/pbeam
      Parameter (idxi = 520)          ! xsi=-ln(xp)
      Parameter (idra = 530)          ! rapidity w.r.t. thrust axis
      Parameter (idpi = 540)          ! Pt in event plane
      Parameter (idpo = 550)          ! Pt out of event plane
      Parameter (offinc = 0)          ! offset

C -- angle of the thrust axis
      Integer idt1,
     >        idt2,
     >        idt3
      Integer offthr

      Parameter (idt1 = 610)          ! thrust angle, before selection
      Parameter (idt2 = 620)          !               after
      Parameter (idt3 = 630)          !               with good chg tracks only
      Parameter (offthr = 0)          ! offset


C -- return values, couting
      Integer idc1,
     >        idc2
      Integer offcou

      Parameter (idc1 = 710)          ! return values
      Parameter (idc2 = 720)          ! couting
      Parameter (offcou = 0)          ! offset


C -- some variables for monitoring
      Integer idet,
     >        idnt,
     >        idfl
      Integer offmcl

      Parameter (idet = 810)          ! total energy of all selected objects
      Parameter (idnt = 820)          ! number of selected objects
      Parameter (idfl = 830)          ! flavour of current event : uds=1,c=2,b=3
      Parameter (offmcl = 0)          ! offset

