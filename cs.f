c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Covid-90 epidemic simulator, april 25 2020
c     S-E1-E2-E3-E4-I-W-B-R where E2,E3 and I are infecterious, 
c                           and B=ICU beds of average of 12 days.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common/sizes1/ worksize,famsize,frsize,crossize,schoolsize
      common/weights/ workinf,school,faminf,frinf,crosinf,ainf(9,9)
      common/weights1/ grandpar
      common/param1/ Nsize,konec,kmax
      common/society/ iage(1500000),iwork(200000,0:50),
     q                ifam(0:1000000,0:12),
     q                ifriend(0:500000,0:20)
      common/society1/ jwork(1500000),
     q            jfam(1500000),jfriend(1500000),jcross(1500000)
      common/society2/ Ncomp,Nschool
      common/links/ Kon(1100000),K(1100000,100)
      common/links1/ Kiden(1100000,100)
      common/strength/ SK(1100000,100)
      common/ages/ jgroup(9,0:200000)
      common/persons/ eta(1500000),agescal(10),actpers(1500000) 
      common/persons1/ actpmax,etamax,agemaxx,actpave,etaave
      common/persons2/ beta0,dt,actu,actp,unspec
      common/persons3/ agescal0(10)
      common/mossong/ polymod(9),polymodg(5)
      common/hospitals/ hospital(10),death(10)
      common/hospitals1/ agecount(10),bedsicu(10)
      common/super1/ superfreq,superpow,cut
      common/super2/ Nsupermax,Nsuperlim,Nsuper(1500000)
      common/super3/ histsuper(0:10000)
      common/super4/ ssuper(1500000)
      common/bbc/ datamat(7,0:10,0:10)
      common/hospitalization/ lockdown1,inhosp(-300:1000)
      common/sweeden/ lockdowns1,insw(-300:1000),ideaths(-300:1000)
      dimension N(1500000)
      dimension tf(-10000:10000),tserial(-10000:10000)
      dimension its(1500000),ipointer(1500000),iserial(1500000)
      dimension ts(-10000:10000)
      dimension itfirst(1500000),istart(1500000)
      dimension gamma(150)
      dimension hist(0:550),histinf(0:550),histw(0:550),histpos(0:550)
      dimension infected(1500000)
      dimension infected0(1500000)
      dimension statage(10),statcross(10,10)
      dimension statsub(8,10,10)
      dimension statgroup(10)
      dimension attack(10),sick(10),hos(10),fat(10)
      dimension resp(10)
      dimension work(0:100)
      dimension isource(1500000)
      dimension nacc(0:10000),iacc(0:10000),nncc(0:10000)
      dimension trajech(-200:3000),trajecicu(-200:3000)
      dimension trajw(-200:3000),trajoth(-200:3000)
      dimension trajsick(-200:3000),trajdead(-200:3000)
      dimension icontraj(-200:3000)
      dimension histacc(0:10000)
      dimension chi2(-50:50)
      dimension ntrajech(-200:3000),ntrajdead(-200:3000)
      dimension atrajech(-200:3000),xtrajech(-200:3000)
      dimension ytrajech(-200:3000)
      dimension xxtrajech(-200:3000),yytrajech(-200:3000)
      dimension atrajdead(-200:3000),xtrajdead(-200:3000)
      dimension ytrajdead(-200:3000)
c
      idum=199729
      call data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     gfortran -fcheck=bounds r.f -o r
c     gfortran -mcmodel=medium -fcheck=bounds r.f -o r
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      chibest=9988888
c
      do i=-200,1000
          ntrajech(i)=0
          ntrajdead(i)=0
          atrajech(i)=0
          xtrajech(i)=0
          xxtrajech(i)=0
          ytrajech(i)=99899.
          yytrajech(i)=99899.
          atrajdead(i)=0
          xtrajdead(i)=0
          ytrajdead(i)=99899.
      enddo
c
      fudgemin=1.
      fudgemax=0
      fudgesmin=1.
      fudgesmax=0
      conmin=1
      conmax=0
      nparam=1
      facacc=5.
c     parameter scan for fitting to data
c     default is that we only use one set of parameters
      do 889 iparam=1,nparam
      do 887 jparam=1,1
      do 886 kparam=1,1
      chisum=0
      chisum2=0
      avekk=0
      kkkmax=0
      if(iparam.gt.1) kkkmax=0
      do 888 kkkc=0,kkkmax
c
      sickfree=0
      sickdays=0
      atest=0.
      iamax=0
      attacmax=0
C     Basic system size and simulation time parameters:
      Nsize=1000000
      size=Nsize
      nday=24
      dt=1./float(nday)    !timestep in unit of days.
      nt=150./dt
c
c.....size of groups to connect to (-1 because self membership)
      worksize=8.-1 
      schoolsize=18.-1
      famsize=2.3-1
      crossize=10.      !normally 10
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c.....Basic supnetworks for infection:
c     unspec=0 correspond to set the random encounters=0
c     workinf=0 correspond to eliminate work contacts
c     feminf=0.2 assign 20% of standard interaction within household
c     crosinf=0.1 means that 10% of interactios go to wider social network
c
c     
c
c.....Default parameters, not to be changed (see later, search for aa)
c     secure that different sectors equally social: 
      unspec=0.82
      faminf=0.54
      grandpar=0
      crosinf=0
      workinf=0.89
      school=0.89   !
c
      unspec0=unspec
      faminf0=faminf
      grandpar0=grandpar
      crosinf0=crosinf
      workinf0=workinf
      school0=school
c
c.....Super spreader option, Nsupermax=inf. activity of superspreader relative to base
c.....                       superfreq=probability that an agent is assigned to be superspreader
c.....                       As Nsupermax is increased, then bet0 should be decreased 
c.....                                                  to fix growth rate of disease.
c.....10% doing 50 times more
      Nsupermax=50
      superfreq=0.10
      Nsupermax=1
      beta0=1.0
      superfreq=0
      superpow=-2.5
c.....Normal parameters
c
      Nsupermax=1
      superfreq=0.
      beta0=1.1
      helpfac=1.
c.....Gamma distributed superspreader
      beta0=1
      Nsuperlim=5000
      superpow=0.9    !=1-k, with k=formfactor in gamma distribution
      Rep0=0.9        !superpow=0.8, Rep0=0.9, cv**2=5.4
      Rep0=1.1        !superpow=0.9, Rep0=1.1
      cut=Rep0/(1.-superpow)
      helpfac=10
c     write(72,*) '#cut=',cut,' if cut=xk what then beta0=',betacorr/Rep0
c.....aaa
c.....
c     option for child less infectious (child=2 means child half infectious)
      child=1
c     option for child less exposed (child=2 means child half exposed)
      childinf=1
c
c     superpow<0 then dont call power laws/gamma distribution
c     0=<superpow<1: gamma distribution
c     superpoe>1 power law
c
c.....initialize social setting and rates
      call initialize(idum)
c.....Call gamma distribution, output Nsuper() and ssuper()
      if(superpow.gt.0) call sspower(idum)
c.....rates adjusted for group sizes, and stored 
      workinf0=workinf    !cluster                 (0.3)
      school0=school
      faminf0=faminf      !cluster                 (0.2)
      grandpar0=grandpar  !cluster to two clusters (0.1)
      crosinf0=crosinf    !networks of friends     (0.1)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c.....to colect statistics on how age groups interact
      do i=1,Nsize
          istart(i)=0
          itfirst(i)=0
          infected(i)=0
          infected0(i)=0
          isource(i)=0
      enddo
      do i=1,10
      statage(i)=0
      do j=1,10
          statcross(i,j)=0
          do kk=1,8
              statsub(kk,i,j)=0
          enddo
      enddo
      statgroup(i)=0
      enddo
c
c.....social network connectivity distribution:
      do i=0,550
         hist(i)=0 
         histw(i)=0 
         histinf(i)=0 
         histpos(i)=0 
      enddo
      do i=1,Nsize
         kk=Kon(i)
         hist(kk)=hist(kk)+1
      enddo
      do i=0,150
        write(8,*) i,hist(i) 
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     BASIC DISEASE PARAMETERS:
c.....distribution of infection rates
c
c.....latency time, IE steps that each are exponential
c     basic disease parameters:
      IE=4
c     N()=0 =susceptible
C     N()=1 =E1     !exposed 1
C     N()=2 =E2     !exposed 2
C     N()=3 =E3     !infectious, presymptomatic
C     N()=4 =E4=IE  !infectious presymptomatic
C     N()=IE+1      !infectious
C     N()=IE+2      !hospital before ICU (not infecting)
C     N()=IE+3      !ICU, say 12 days
C     N()=IE+4      !removed
C
      atime=0
      do i=1,IE
          gamma(i)=0.8            !rate of leaving each of IE presymptomatic state
          atime=1/gamma(i)+atime
      enddo
      agamma=1./atime
      delta=0.3333     !rate to leave I state (IE=5)
      hosp=0.5         !time to hospitalization after infecterous
cwww here we made a delay from symptom to hospitalization, now 9 days instead of 5
      hospstay=0.2     !time to hospitalization
      serious=1./12.   !12 days in ICU in Norway
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c.....start condition with 10 infected people (age between 21 and 40)
      do i=1,Nsize
         N(i)=0
      enddo
      Nsize9=Nsize/9
      do i1=1,100   !50  !50
         i=2*Nsize9+4*ran2(idum)*Nsize9
         N(i)=IE+1
      enddo
c
c.....Quarantine Strategies:
c
      icontact=0          !=1, only person, =2 also persons connections
      contact=0.0         !0.1 or delta     !rate of contact tracing
      intcheck=200./dt
      ncheck=0 
c.....
c
c.....write out parameters:
      R0=abeta/delta
      alow=0
      ahigh=0
      write(6,110) abeta,dbeta,agamma,delta,IE,R0,
     q             alow,ahigh,icontact,contact,
     q             intcheck,ncheck,ak
c     write(16,110) abeta,dbeta,agamma,delta,IE,R0,
c    q             alow,ahigh,icontact,contact,
c    q             intcheck,ncheck,ak
c     write(26,110) abeta,dbeta,agamma,delta,IE,R0,
c    q             alow,ahigh,icontact,contact,
c    q             intcheck,ncheck,ak
 110  format('# beta,gamma,delta=',4F7.3,I4,F7.2,
     q       ' pop=',2F6.3,' icontact=',I2,F6.3,' test=',2I4,
     q       ' ave k=',3F7.2)
c.....
      do i=0,10000
         tf(i)=0
         tf(-i)=0
         tserial(i)=0
         tserial(-i)=0
      enddo
c.....
c
c.....special analysis for super spreaders:
      itest=0
      testlim=0.005
      do i=0,100 
         work(i)=0
      enddo
c
c.....basic analysis
      sumspec=0
      sumunspec=0
      sttac=0
      hospinl=0
      bedsmax=0
      sickmax=0
      hospmax=0
      testmax=0
      bedssum=0
      infsum0=0
      xsum0=0
      fatsum=0
      superinf=0
      suminf=0
      suminf2=0
      baseinf=0
c
      countf=0
      earlyinf=0
      infectf=0
      infectf2=0
      noinfectf=0
      already=0
      already2=0
      do i=0,10000
         iacc(i)=0
         nacc(i)=0
         nncc(i)=0
      enddo
      do i=-200,3000
          trajoth(i)=1
          trajw(i)=1
          icontraj(i)=-1
          trajech(i)=0
          trajecicu(i)=0
          trajsick(i)=0
          trajdead(i)=0
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     actual parameters to run with:
czzzaaa
      unspec=unspec0
      faminf=faminf0      !cluster                 (0.2)
      workinf=workinf0    !cluster                 (0.3)
      school=school0
      grandpar=grandpar0  !cluster to two clusters (0.1)
      crosinf=crosinf0    !networks of friends     (0.1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c.....specific has weight spec*actpave, spec=friends if no Mitigation
c.....public has weight   unspec*etaave, unspec=1-friends if no M.
      actu=unspec*etaave
      acts=actpave
      write(72,169) actu,acts,
     q              etamax,actpmax,etaave,actpave
      write(41,169) actu,acts,
     q              etamax,actpmax,etaave,actpave
      write(71,169) actu,acts,
     q              etamax,actpmax,etaave,actpave
      actu=actu*dt
      acts=acts*dt
      write(72,169) actu,acts,
     q              actu*helpfac,acts*helpfac
 169  format('# actu,acts=',2F7.3,' etam,actp=',4F8.3)
      acturemove=0
      actunorm=0
      actsremove=0
      actsnorm=0
c
c.....if containlim1>1 then no mitigation is attempted
c.....otherwise mitigitation start when accumulated cases containlim1
      icon=0
      containlim1=0.01
      tcountlim1=30
      tcountlim2=70
      tcountlim3=200
      tcountlim1=200
      tcount=0
      detailmax=0
c.....kstat=0 is standard, kstat= is only for samling group activity
c.....independent of disease parameters
      fudge0=1.  
      fudges0=1.0
      kstat=0
      xaaa=0
      yaaa=0
      ymean=0
      yvar=0
      do 999 it=1,nt
c.........to measure R0, R:
          if(nday*(it/nday).eq.it) then
              infsum=0
              xsum=0
          endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         progress of disease for all (infected) agents
c.........id=1,2,3,...IE,IE+1 (controlled by gamma)
c
          hospinl=0
          do i=1,Nsize
             id=N(i)
             if(id.ge.1.and.id.le.IE) then
                 if(ran2(idum).lt.gamma(id)*dt) N(i)=id+1
                 if(N(i).eq.IE+1) then
                     iii=it-itfirst(i)
                     tf(iii)=tf(iii)+1
                     istart(i)=it
                 endif
             endif
          enddo
c
c.........Disease can be cured
c.........id=IE+1,IE+2 (controlled by delta)
c.........id=IE+2,IE+3 (controlled by hosp)
c.........id=IE+3,IE+4 (controlled by serious=1/12 in Noway data)
c
          do ig=1,10
             hos(ig)=0    !hospital arrival in age group ig
             fat(ig)=0    !fatality in age group ig
          enddo
          do i=1,Nsize
             id=N(i)
             ig=iage(i)
             if(id.eq.IE+1) then
                 if(ran2(idum).lt.delta*dt) then
                     N(i)=IE+2
                     xaaa=xaaa+1
                     xsum=xsum+1
                     xsum0=xsum0+1
                     infsum=infsum+infected(i)
                     iiii=infected(i)
                     jjjj=Nsuper(i)
                     histacc(jjjj)=histacc(jjjj)+iiii
                     if(sttac.gt.0.05*Nsize.and.
     q                  sttac.lt.0.09*Nsize) then
                         yaaa=yaaa+1
                         ymean=ymean+iiii
                         yvar=yvar+iiii**2
                         nncc(iiii)=nncc(iiii)+1
                         iacc(iiii)=iacc(iiii)+iiii
                     endif
                     infsum0=infsum0+infected0(i)
                     infsum00=infsum00+infected0(i)**2
                 endif
             else if(id.eq.IE+2) then
                 if(ran2(idum).lt.hosp*dt) then
                 N(i)=IE+3
                 hos(ig)=hos(ig)+hospital(ig)
                 hospinl=hospinl+hospital(ig)*agecount(ig)
                 endif
             else if(id.eq.IE+3) then
                 if(ran2(idum).lt.hospstay*dt) then
                 N(i)=IE+4
                 ig=iage(i)
                 fat(ig)=fat(ig)+1.0*death(ig)
cpeople below 70 go t0 ICU
cpeople above 70 die. i.e. approximate with deth bedor ICU 
                 endif
             else if(id.eq.IE+4) then
                 if(ran2(idum).lt.serious*dt) then
                     N(i)=IE+5
                     ig=iage(i)
                     fat(ig)=fat(ig)+0.*death(ig)
                 endif
             endif
          enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Main infection loop, done with 
c         propability=activity for each agent: 
c
          do 998 index=1,Nsize
c
c......   Agent infect specifically:
c......   do infections with probability = acts 
c......   (and then chose person with prop actpers(i))
c
          if(ran2(idum).lt.acts*helpfac) then
c......          chose a person to infect from
 83              i=ran2(idum)*Nsize+1
                 ig=iage(i)
                 if(actpers(i).gt.1) write(6,*) 'wwwp'
c................actpers() is normalized to maximum, just chose 
c................one agent here weighted with actpers()
                 if(ran2(idum).gt.actpers(i)) goto 83
                 id=N(i)
                 indc=0
c......          Need to be infectious
                 if(id.eq.IE+1.or.id.eq.IE.or.
     q              id.eq.IE-1.or.kstat.eq.1) then
c.................................................................
c
c                loop in case agen i is a super spreader:
c
                 do 990 isuper=1,Nsuper(i)
                 if(ran2(idum)*Nsuper(i)*helpfac.lt.ssuper(i)) then  
c                in case of child>1 the children (ig=1,ig=2) is less infectious                         
                 if(ig.gt.2.or.child*ran2(idum).lt.1) then
 85              j1=ran2(idum)*Kon(i)+1
                 j=K(i,j1)
                 if(ran2(idum).gt.SK(i,j1)) goto 85
                 jg=iage(j)
                 Kii=Kiden(i,j1)
c
c................Age profile, less active persons will be less infected
                 actsnorm=actsnorm+1
                 if(ran2(idum).gt.agescal(jg)) goto 991
                 actsremove=actsremove+1
c
c................Heterogenous agents, the minimum sets the interaction:
                 if(eta(j).lt.eta(i)) then
                     zzz=eta(j)/eta(i)
                     if(ran2(idum).gt.zzz) goto 991
                 endif
c................Disease state of agent j:
                 jd=N(j)
                 jg=iage(j)
                 if(jg.le.2.and.childinf*ran2(idum).gt.1) goto 991
                     if(Kiden(i,j1).eq.1) then
                        if(ran2(idum)*faminf0.gt.faminf) goto 991
                     endif
                     if(Kiden(i,j1).eq.2) then
                        if(ran2(idum)*grandpar0.gt.grandpar) goto 991
                     endif
                     if(Kiden(i,j1).eq.3) then
                        if(ran2(idum)*crosinf0.gt.crosinf) goto 991
                     endif
                     if(Kiden(i,j1).eq.4) then
                        if(ran2(idum)*workinf0.gt.workinf) goto 991
                     endif
                     if(Kiden(i,j1).eq.5) then
                        if(ran2(idum)*school0.gt.school) goto 991
                     endif
c................Sampling of Matrix of subgroup interactins:
                     if(Kii.eq.1) then
                        statgroup(1)=statgroup(1)+1
                        statsub(1,ig,jg)=statsub(1,ig,jg)+1
                     else if(Kii.eq.2) then
                        statgroup(2)=statgroup(2)+1
                        statsub(2,ig,jg)=statsub(2,ig,jg)+1
                     else if(Kii.eq.3) then
                        statgroup(3)=statgroup(3)+1
                        statsub(3,ig,jg)=statsub(3,ig,jg)+1
                     else if(Kii.eq.4) then
                        statgroup(4)=statgroup(4)+1
                        statsub(4,ig,jg)=statsub(4,ig,jg)+1
                     else if(Kii.eq.5) then
                        statgroup(5)=statgroup(5)+1 
                        statsub(5,ig,jg)=statsub(5,ig,jg)+1
                     endif
c
c....................lock down/strategu simulation in special periods
c....................
                     sumspec=sumspec+1
                     statage(ig)=statage(ig)+1
                     statage(jg)=statage(jg)+1
                     statcross(ig,jg)=statcross(ig,jg)+1
                     if(isource(j).ne.i) then
                     infected0(i)=infected0(i)+1      !to measure R0
                     endif
c
                     if(jd.gt.0) then
                         if(attac.lt.0.05) then
                           if(Kii.eq.1) then
                              already=already+1
                              if(ifam(inum,0).eq.2) then
                                 already2=already2+1
                              endif
                           endif
                         endif
                     endif
                     if(jd.eq.0) then
                         infected(i)=infected(i)+1    !to measure R
                         isource(j)=i
                         N(j)=1
                         itfirst(j)=it
                         isss=Nsuper(i)
c                        time from when i infected to subsequent infection
                         its(j)=it-itfirst(i)
                         ipointer(j)=i
                         sttac=sttac+1
                         x0sum=x0sum+1
                         xksum=xksum+Kon(j)
                         suminf=suminf+1
                         if(attac.lt.0.05) then
                           earlyinf=earlyinf+1
                           inum=jfam(i) 
                           if(ifam(inum,0).eq.2) then
                                suminf2=suminf2+1
                           endif
                           if(Kii.eq.1) then
                              infectf=infectf+1
                              if(ifam(inum,0).eq.2) then
                                   infectf2=infectf2+1
                              endif
                           else
                              noinfectf=noinfectf+1
                           endif
                         endif
                         if(Nsuper(i).gt.1) superinf=superinf+1
                         if(Nsuper(i).eq.1) baseinf=baseinf+1
                     endif
c                this was end of infection 
 991             continue
                 endif
                 endif
 990             continue
                 endif
             endif
c
c......      Now agent ingfect publically:
c......      Do infections with probability=actu
             if(ran2(idum).lt.helpfac*actu) then
             if(actu.gt.1) write(6,*) 'wwwu'
c......          chose person to infect from (with prop eta*agescal)
 84              i=ran2(idum)*Nsize+1
                 ig=iage(i)
c......          
                 if(eta(i)*agescal(ig)/etamax.gt.1) write(6,*) 'uu'
                 if(ran2(idum)*etamax.gt.
     q              eta(i)*agescal(ig)) goto 84
                 id=N(i)
                 indc=0
c......          Need to be infecterous
                 if(id.eq.IE+1.or.id.eq.IE.or.
     q              id.eq.IE-1.or.kstat.eq.1) then
                 if(indc.eq.0) then
                 do 993 isuper=1,Nsuper(i)
                 if(ran2(idum)*Nsuper(i)*helpfac.lt.ssuper(i)) then  
c................in case we have a child it may be less infectious 
c                (with factor ``child", when child>1):
                 if(ig.gt.2.or.child*ran2(idum).lt.1) then
c....................public infections (with age and quenched activity):
 86                  j=ran2(idum)*Nsize+1
                     jg=iage(j)
                     if(jg.le.2.and.childinf*ran2(idum).gt.1) goto 992
                     actunorm=actunorm+1
                     if(ran2(idum).gt.agescal(jg)) goto 992
                     if(eta(j).lt.eta(i)) then
                         zzz=eta(j)/eta(i)
                         if(ran2(idum).gt.zzz) goto 992
                     endif
                     acturemove=acturemove+1
                     sumunspec=sumunspec+1
                     statage(ig)=statage(ig)+1
                     statage(jg)=statage(jg)+1
                     statcross(ig,jg)=statcross(ig,jg)+1
                     statgroup(6)=statgroup(6)+1
                     statsub(6,ig,jg)=statsub(6,ig,jg)+1
                     jd=N(j)
                     if(isource(j).ne.i) infected0(i)=infected0(i)+1
                     if(jd.eq.0.or.kstat.eq.1) then
                         infected(i)=infected(i)+1
                         isource(j)=i
                         N(j)=1
                         isss=Nsuper(i)
                         nacc(isss)=nacc(isss)+1
                         itfirst(j)=it
                         its(j)=it-itfirst(i)
                         ipointer(j)=i
                         sttac=sttac+1
                         x0sum=x0sum+1
                         xksum=xksum+Kon(j)
                         suminf=suminf+1
                         if(attac.lt.0.05) then
                             earlyinf=earlyinf+1
                             inum=jfam(i) 
                             if(ifam(inum,0).eq.2) then
                                 suminf2=suminf2+1
                             endif
                             noinfectf=noinfectf+1
                         endif
                         if(Nsuper(i).gt.1) superinf=superinf+1
                         if(Nsuper(i).eq.1) baseinf=baseinf+1
                     endif
 992             continue
                 endif
                 endif
 993             continue
                 endif
                 endif
c                end of quarentine 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 995      continue
          endif
 998      continue
c
c         infections during timestep dt finished 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c.........statistics at time it*dt:
          infect=0
          infecf=1
          attac=0
          expos=0
          karan=1
          do ig=1,10
               resp(ig)=0
               sick(ig)=0
               attack(ig)=0
          enddo
          do i=1,Nsize
              ig=iage(i)
              if(N(i).gt.0.and.N(i).le.IE) expos=expos+1
              if(N(i).eq.IE+1.or.N(i).eq.IE+2) then
                  infect=infect+1
                  sick(ig)=sick(ig)+1
              endif
              if(N(i).eq.IE+4) then
c                 people that goes into ICU (stay there 12 days)
                  resp(ig)=resp(ig)+1
              endif
              if(N(i).eq.IE+1) then
                  infecf=infecf+1
              endif
              if(N(i).ge.IE+1) then
                  attac=attac+1
                  attack(ig)=attack(ig)+1
              endif
          enddo
          sickfree=sickfree+infecf*dt
          sickdays=sickdays+infect*dt
          if(infect.gt.attacmax) then
              iamax=it
              attacmax=infect
          endif
          xksum=0
          if(x0sum.gt.0) xksum=xksum/x0sum
          Rnow=0
          if(xsum.gt.0) Rnow=infsum/xsum
          if(xsum0.gt.0) then
             Rnow0=infsum0/xsum0
             Vnow0=0.0 
             if(Rnow0.gt.0) then
             Vnow0=(infsum00/xsum0-Rnow0**2)/Rnow0**2
             endif
          endif
c
c.........fraction of involved population in each age group:
          Nsize9=Nsize/9
          do ig=1,10
              resp(ig)=resp(ig)/Nsize9
              sick(ig)=sick(ig)/Nsize9
              attack(ig)=attack(ig)/Nsize9
              hos(ig)=hos(ig)/Nsize9
              fat(ig)=fat(ig)/Nsize9
          enddo
c
c.........Estimate of young sick people:
          sickyoung=sick(1)*agecount(1)
     q              +sick(2)*agecount(2)
     q              +sick(3)*agecount(3)
          attyoung=attack(1)*agecount(1)
     q             +attack(2)*agecount(2)
     q             +attack(3)*agecount(3)
          attyoung=attyoung/(agecount(1)+agecount(2)+agecount(3))
c
c.........ICU weights from official norwegian model:
          beds=0
          sickall=0
          hospinl=0
          fatality=0
          do ig=1,9
              sickall=sickall+sick(ig)*agecount(ig)       
              beds=beds+resp(ig)*bedsicu(ig)*agecount(ig)  
              hospinl=hospinl+hos(ig)*agecount(ig)
              fatality=fatality+fat(ig)*agecount(ig)
          enddo
          fatsum=fatsum+fatality
          hospinl=hospinl/dt
          fatality=fatality/dt
          fatold=fatality
          hosold=hospinl
          if(beds.gt.bedsmax) bedsmax=beds 
          if(sickall.gt.sickmax) sickmax=sickall 
          if(hospmax.lt.hospinl) hospmax=hospinl 
          bedssum=bedssum+beds*dt
          attac=attac/Nsize
          sickold=sick(8)*agecount(8)+agecount(9)*sick(9)
          sickold=sickold/(agecount(8)+agecount(9))
          attold=attack(8)*agecount(8)+attack(9)*agecount(9)
          attold=attold/(agecount(8)+agecount(9))
          if(iparam.eq.1.and.kkkc.le.3) then
          trajwww=(5*workinf+2*school)/(5*workinf0+2*school0)
          trajooo=unspec/unspec0
          write(6,125) it*dt,fatality,sickall,
     q               attac,sickyoung,sick(7),
     q               sickold,
     q               sickall,hospinl,beds,attyoung,
     q               trajwww,trajwww+trajooo,attold,
     q               Rnow,Rnow0,icon,Vnow0
          itt=it*dt
          trajw(itt)=(5*workinf+2*school)/(5*workinf0+2*school0)
          trajoth(itt)=unspec/unspec0
          endif
          itt=it*dt
          icontraj(itt)=icon
          pop=5728769
          trajech(itt)=trajech(itt)+hospinl*pop*dt
          trajecicu(itt)=trajecicu(itt)+beds*pop*dt
          trajsick(itt)=trajsick(itt)+sickall*pop*dt
          trajdead(itt)=trajdead(itt)+fatality*pop*dt
c
c.........special testing for superspreaders
          if(attac.gt.testlim.and.itest.eq.0) then
              itest=1
              sumh=0
              do inum=Nschool+1,Nschool+Ncomp
                 infec=0
                 ipos=0
                 isum=0 
c                if(iwork(inum,0).gt.5.and.iwork(inum,0).le.15) then
                 if(iwork(inum,0).gt.0) then
                 do kk=1,iwork(inum,0)
                    i=iwork(inum,kk)
                    ig=iage(i)
                    id=N(i)
                    if(id.ge.IE-1.and.id.le.IE+2) infec=infec+1
                    if(id.ge.IE-1) ipos=ipos+1
                    isum=isum+1
                 enddo
                 endif
                 histw(isum)=histw(isum)+1
                 histinf(infec)=histinf(infec)+1
                 histpos(ipos)=histpos(ipos)+1
                 sumh=sumh+1
              enddo 
              do ix=0,150
                 histw(ix)=histw(ix)/sumh
                 histinf(ix)=histinf(ix)/sumh
                 histpos(ix)=histpos(ix)/sumh
              enddo
          endif
c.........end statistics
c
c.........Oossibility to impose special rules
c         when attac is large than som limit, and some 
c         time passed after that (here 30 days)
c
          if(attac.gt.containlim1.and.icon.eq.0) then
c             3.5 million-->1.6million->2.1million work
c             1.3 million-->0         ->0.5million
c             4.8 million-->1.6million->2.6million
              write(72,372) it*dt,attac,sickall 
              modellock1=it*dt
 372          format('# t=',F8.1,' attac,sick=',2F9.4)
              fudge=fudge0
              fudges=fudges0
              icon=1
              tcount=0
              unspec=unspec0*fudges
              faminf=faminf0
              workinf=workinf0*fudge    !real data, 0.45
              school=school0*fudge      !real data, 0.0
              grandpar=grandpar0*fudge
              crosinf=crosinf0*fudge
              actu=unspec*etaave*dt
              acts=actpave*dt
              xsum0=0
              infsum0=0
              do i=1,9
                 agescal(i)=agescal0(i)
              enddo
c             agescal(8)=agescal0(8)/10
c             agescal(9)=agescal0(9)/10
          endif
          if(tcount.gt.tcountlim1.and.icon.eq.1) then
              write(72,372) it*dt,attac,sickall 
              fudge=fudge0
              fudges=fudges0
              icon=2
              unspec=unspec0*fudges
              faminf=faminf0
              workinf=workinf0*fudge   !*0.60  !real data, physical work fraction in april
              school=school0*fudge    !*0.45
              grandpar=grandpar0*fudge  
              crosinf=crosinf0*fudge
              actu=unspec*etaave*dt
              acts=actpave*dt
              xsum0=0
              infsum0=0
              do i=1,9
                 agescal(i)=agescal0(i)
              enddo
c             agescal(8)=agescal0(8)/10
c             agescal(9)=agescal0(9)/10
          endif
          if(tcount.gt.tcountlim2.and.icon.eq.2) then
              write(72,372) it*dt,attac,sickall 
              fudge=fudge0
              fudges=fudges0
              icon=3
              unspec=unspec0*fudges
              faminf=faminf0
              workinf=workinf0*fudge
              school=school0*fudge
              grandpar=grandpar0*fudge
              crosinf=crosinf0*fudge
              actu=unspec*etaave*dt
              acts=actpave*dt
              xsum0=0
              infsum0=0
              do i=1,9
                 agescal(i)=agescal0(i)
              enddo
          endif
          if(tcount.gt.tcountlim3.and.icon.eq.3) then
              write(72,372) it*dt,attac,sickall 
              fudge=fudge0   
              fudges=1
              icon=4
              unspec=unspec0*fudges
              faminf=faminf0
              workinf=workinf0*fudge
              school=school0*fudge
              grandpar=grandpar0*fudge
              crosinf=crosinf0*fudge
              actu=unspec*etaave*dt
              acts=actpave*dt
              xsum0=0
              infsum0=0
              do i=1,9
                 agescal(i)=agescal0(i)
              enddo
          endif
          tcount=tcount+dt
 999  continue
 125      format(' ',F7.2,F10.7,6F9.6,2F12.8,F11.5,3F8.5,
     q               2F6.2,I2,2F8.4)
 126      format('# Attac fraction=',F9.4,
     q           ' max number in IE+1,IE+2=',F9.4,
     q           ' # max ICU beds=',F10.6,
     q           ' total ICU bed days=',F9.6,
     q           ' indl=',F9.6,' Death=',2F9.6)
      acturemove=acturemove/actunorm
      actsremove=actsremove/actsnorm
      weightx=(actu*acturemove+acts*actsremove)
      weightr=5.5*(actu*acturemove+acts*actsremove)/dt
      weights=(actu*acturemove+acts*actsremove)/(acts+actu)
      timeinf=1/delta+1/gamma(IE)+1/gamma(IE-1)
      write(72,168) acturemove,actsremove,weightx,
     q              weightr,weightr/beta0,weights,
     q              timeinf,dt
 168  format('# removal factors: actu,acts=',2F7.3,
     q       ' weights=',2F8.3,' weight/beta0=',
     q         F8.3,' chance to remove=',F8.3,' tinf=',2F6.2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     distribution of latency times:
      if(kkkc.eq.0) then
      xserial=0
      do j=1,Nsize
          i=ipointer(j)
          if(i.gt.0) then
          jdif=istart(j)-istart(i)
          tserial(jdif)=tserial(jdif)+1
          xserial=xserial+1
          endif
      enddo
      xtime=0.1
      atime=0
      stime=0.1
      xstime=0
      ctime=0.1
      astime=0
      do i=-10000,10000
          ts(i)=0
      enddo
      sum=0.01
      sums=0.01
      do i=1,Nsize
          if(its(i).gt.0) then
             sum=sum+1
             ii=its(i)
             ts(ii)=ts(ii)+1
          endif
      enddo
      do i=-1000,1000
          if(ts(i).gt.0) then 
              xstime=xtime+ts(i)
              astime=astime+i*ts(i)*dt
          endif
          if(tf(i).gt.0) then 
              xtime=xtime+tf(i)
              atime=atime+i*tf(i)*dt
          endif
          if(tserial(i).gt.0) then 
              sums=sums+1
              ctime=ctime+tserial(i)
              stime=stime+i*tserial(i)*dt
          endif
          if(tf(i).gt.0.or.ts(i).gt.0.or.tserial(i).gt.0) then
          write(7,*) i*dt,tf(i)/(attac*dt),
     q              ts(i)/(sum*dt),tserial(i)/(xserial*dt)
          endif
      enddo
c.....average latency time (should=5), average serial time (7.5):
      write(7,120) atime/xtime,astime/xstime,stime/ctime
      endif
      sickfree=sickfree/size
      sickdays=sickdays/size
      carandays=0.
      karanday=0
      testmax=testmax/dt
      if(kkkc.ge.0) then
      write(12,130) icontact,karanday,contact,attac/Nsize,
     q             attacmax/size,int(iamax*dt),
     q             sickdays,carandays,sickfree,atest,
     q             ncheck,intcheck,bedsmax,hospmax,testmax
      else
      write(13,130) icontact,karanday,contact,attac/Nsize,
     q              attacmax/size,int(iamax*dt),
     q              sickdays,carandays,sickfree,atest,
     q              ncheck,intcheck,bedsmax,hospmax,testmax
      endif
c
      do i=1,550
          hist(i)=0
      enddo
      sum=0
      do i=1,Nsize
         iii=infected(i)
         ig=iage(i)
         sum=sum+agecount(ig)
         if(iii.le.550) hist(iii)=hist(iii)+agecount(ig)
      enddo
      prop=1
      do i=1,550
         x=i
         if(hist(i)+histinf(i)+histw(i).gt.0) then
         write(81,181) i,hist(i)/sum,
     q                 histinf(i),histw(i),histpos(i)
         endif
      enddo
 181  format(' ',I5,6F10.6)
c
      if(iparam.eq.1) then
      idiff=modellock1-lockdown1
      idiffs=modellock1-lockdowns1
      endif
      idiff=idiffs
      idiff=0
      write(78,301) lockdown1,modellock1,idiff,
     q            lockdowns1,modellock1,idiffs
      write(180,301) lockdown1,modellock1,idiff,
     q            lockdowns1,modellock1,idiffs
      ntt=nt*dt
c     ratio of sweedish to danish population
      ratpop=5790037./10092558.
      do kc=-10,30
          chi2(kc)=0
      enddo
      do i=-3,ntt-idiff
c
      write(180,881) iparam,kkkc,Nsupermax,superfreq,
     q               fudge0,fudges0,containlim1 
 881  format('# parameters=',3I5,F6.2,' fudges=',3F9.5)
c.....
      ratsw=5790037./10092558.
      if(ideaths(i).ge.0) then
      do kc=-10,20
          ikc=i+idiff+kc
          expected=trajdead(i+idiff+kc)/ratsw
          xxx=0
          if(expected.gt.5) then
              xxx=(ideaths(i)-expected)**2/expected
          endif
          chi2(kc)=chi2(kc)+xxx
      enddo
      endif
      if(iparam.le.1) then
      if(trajech(i+idiff).gt.xxtrajech(i+idiff)) then
          xxtrajech(i+idiff)=trajech(i+idiff)
      endif
      if(trajech(i+idiff).lt.yytrajech(i+idiff)) then
          yytrajech(i+idiff)=trajech(i+idiff)
      endif
c.....
      if(ideaths(i).ge.0) then
      if(kkkc.eq.0) then
      write(78,878) i,inhosp(i),trajech(i+idiff),
     q            trajech(i),trajecicu(i+idiff),
     q            trajw(i+idiff),
     q            trajw(i+idiff)+trajoth(i+idiff),
     q            insw(i+idiff),trajsick(i+idiff),
     q            trajdead(i+idiff),ideaths(i)
      else
      write(79,878) i,inhosp(i),trajech(i+idiff),
     q            trajech(i),trajecicu(i+idiff),
     q            trajw(i+idiff),
     q            trajw(i+idiff)+trajoth(i+idiff),
     q            -999,trajsick(i+idiff),
     q            trajdead(i+idiff),-999
      endif
      else
      if(kkkc.eq.0) then
      write(78,878) i,inhosp(i),trajech(i+idiff),
     q            trajech(i),trajecicu(i+idiff),
     q            trajw(i+idiff),
     q            trajw(i+idiff)+trajoth(i+idiff),
     q            -999,trajsick(i+idiff),
     q            trajdead(i+idiff)
      endif
      if(kkkc.gt.0) then
      write(79,878) i,inhosp(i),trajech(i+idiff),
     q            trajech(i),trajecicu(i+idiff),
     q            trajw(i+idiff),
     q            trajw(i+idiff)+trajoth(i+idiff),
     q            -999,trajsick(i+idiff),
     q            trajdead(i+idiff),-999
      endif
      endif
c.....
      endif
c
      enddo
c
      chimin=999999
      do kc=-10,20
      if(chi2(kc).lt.chimin) then
           chimin=chi2(kc)
           kkmin=kc
      endif
      enddo
      chisum=chisum+chimin
      chisum2=chisum2+chimin**2
      avekk=avekk+kkmin
      write(78,879) kkmin,chimin,0.,fudge0,
     q              fudges0,containlim1,iparam
      write(978,879) kkmin,chimin,0.,fudge0,
     q              fudges0,containlim1,iparam,kkkc
c
      if(iparam.eq.1) then
c     if(chimin.lt.facacc*chibest) then
      write(978,893) chibest,chimin,fugge0,fudges0,containlim1
 893  format('# included into trajectories',2F7.1,
     q       ' par=',3F10.4)
      do i=-3,ntt-idiff
      ntrajech(i+idiff)=ntrajech(i+idiff)+1
      atrajech(i+idiff)=atrajech(i+idiff)+trajech(i+idiff)
      if(trajech(i+idiff).gt.xtrajech(i+idiff)) then
          xtrajech(i+idiff)=trajech(i+idiff)
      endif
      if(trajech(i+idiff).lt.ytrajech(i+idiff)) then
          ytrajech(i+idiff)=trajech(i+idiff)
      endif
c     if(iparam.eq.2) ytrajech(i+idiff)=trajech(i+idiff)
c     if(iparam.eq.3) xtrajech(i+idiff)=trajech(i+idiff)
c
      ntrajdead(i+idiff)=ntrajdead(i+idiff)+1
      atrajdead(i+idiff)=atrajdead(i+idiff)
     q                   +trajdead(i+idiff)
      if(trajdead(i+idiff).gt.xtrajdead(i+idiff)) then
          xtrajdead(i+idiff)=trajdead(i+idiff)
      endif
      if(trajdead(i+idiff).lt.ytrajdead(i+idiff)) then
          ytrajdead(i+idiff)=trajdead(i+idiff)
      endif
      enddo
c     endif
      endif
c
      write(78,*)
      write(79,*)
      do i=-3,ntt-idiff
          icon=icontraj(i+idiff)
          if(icon.eq.1) then
             write(178,*) i,0
          else if(icon.eq.2) then
             write(278,*) i,0
          else if(icon.eq.3) then
             write(378,*) i,0
          else if(icon.eq.4) then
             write(478,*) i,0
          endif
       enddo
      write(6,*)
 301  format('# lock down real,model, idiff=',3I6,
     q       ' sweden, real, model, idiffs=',3I6)
c
 888  continue
      chisum=chisum/kkkc
      chisum2=chisum2/kkkc-chisum**2
      chisum2=sqrt(chisum2)
      avekk=avekk/kkkc
      if(chisum.lt.chibest) then
          chibest=chisum
          chi2best=chisum2
          f0best=fudge0
          fsbest=fudges
          cbest=containlim1
          bestk=avekk
      endif
      write(78,899) avekk,chisum,chisum2,
     q              fudge0,fudges0,containlim1,
     q              iparam,kkkc
      if(chisum.lt.facacc*chibest) then
      if(fudge0.lt.fudgemin) fudgemin=fudge0
      if(fudge0.gt.fudgemax) fudgemax=fudge0
      if(fudges0.lt.fudgesmin) fudgesmin=fudges0
      if(fudges0.gt.fudgesmax) fudgesmax=fudges0
      if(containlim1.lt.conmin) conmin=containlim1
      if(containlim1.gt.conmax) conmax=containlim1
      write(978,899) avekk,chisum,chisum2,
     q              fudge0,fudges0,containlim1
      endif
 886  continue
 887  continue
 889  continue
      write(978,*) '####result=',nparam,kkc,
     q              superpow,superfreq
      write(78,899) bestk,chibest,chi2best,
     q              f0best,fsbest,cbest
      write(978,899) bestk,chibest,chibest2,
     q              f0best,fsbest,cbest
 878  format(' ',2I6,5F12.3,I7,2F12.1,I7)
 879  format('# chi2=',I4,2F10.1,' fudges=',3F7.4,3I4)
 899  format('# chi2=',F4.1,2F10.1,' fudges=',3F7.4,2I4)
c
      write(180,892) fudgemin,fudgemax,fudgesmin,
     q               fudgesmax,conmin,conmax
      write(978,892) fudgemin,fudgemax,fudgesmin,
     q               fudgesmax,conmin,conmax
 892  format('# fudge0 limits=',4F8.3,' con limits=',2F8.3)
      if(fudge0.gt.fudgemax) fudgemax=fudge
      if(containlim1.lt.conmin) conmin=containlim1
      if(containlim1.gt.conmax) conmax=containlim1
      do i=-3,ntt-idiff
      if(ntrajech(i+idiff).gt.0) then
          atrajech(i+idiff)=atrajech(i+idiff)/ntrajech(i+idiff)
      endif
      if(ntrajdead(i+idiff).gt.0) then
      atrajdead(i+idiff)=atrajdead(i+idiff)/ntrajdead(i+idiff)
      endif
      if(xtrajdead(i+idiff).gt.999) xtrajdead(i+idiff)=0
      if(ytrajdead(i+idiff).le.0) ytrajdead(i+idiff)=0
      if(xtrajech(i+idiff).gt.999) xtrajech(i+idiff)=0
      if(ytrajech(i+idiff).le.0) ytrajech(i+idiff)=0
      dummy=atrajdead(i+idiff)
      write(80,878) i,inhosp(i),atrajech(i+idiff)
c    q            ytrajech((i+idiff),
c    q            xtrajech(i+idiff)
c    q            ytrajdead(i+idiff),
c    q            xtrajdead(i+idiff),
c    q            insw(i+idiff)
c    q            trajsick(i+idiff)
c    q            dummy,
c    q            ideaths(i)
      write(180,872) i,inhosp(i),atrajech(i+idiff),
     q            ytrajech(i+idiff),
     q            xtrajech(i+idiff),
     q            ytrajdead(i+idiff),
     q            xtrajdead(i+idiff),
     q            insw(i+idiff),
     q            xxtrajech(i+idiff),
     q            atrajdead(i+idiff),
     q            ideaths(i)
      enddo
 872  format(' ',2I6,5F12.3,I7,2F12.1,I7)
c
c.....tests:
c
      sum=0.1
      sumc=0.1
      do i=1,9
          sum=sum+statage(i)
          do j=1,9
              sumc=sumc+statcross(i,j)
          enddo
      enddo
c
      write(71,170) sumunspec,sumspec 
      write(71,126) bedsmax,bedssum,hospmax
      write(6,126) attac,sickmax,bedsmax,bedssum,hospmax,
     q             fatsum,fatsum/attac
c
c.....output for partition into various work spheres:
c
      if(suminf.ge.0) then
      rf=0
      if(baseinf.gt.0) rf=superinf/baseinf
      write(88,188) baseinf/suminf,superinf/suminf,rf,
     q              Nsupermax,superfreq,Rnow0
      write(72,188) baseinf/suminf,superinf/suminf,rf,
     q              Nsupermax,superfreq,Rnow0
      endif
 188  format('#basal inf=',3F10.4,' Nsuper=',I4,2F10.3) 
      do i=1,9
      write(88,*) i,attack(i),statage(i)/sum
      enddo
      sumx=0
      do i=1,7
      sumx=sumx+datamat(i,0,0)
      enddo
      write(72,177) 
      write(72,179) 
      write(72,174) Rnow,faminf0,grandpar0,
     q              crosinf0,workinf0,school0,unspec0,school0+workinf0
      write(72,174) Rnow0,faminf0*famsize,
     q              grandpar0*(2*famsize+1),
     q              crosinf0*crossize,workinf0*worksize,school0*worksize*3,
     q              unspec0,school0*worksize*3+workinf0*worksize
      if(earlyinf.gt.0) then
           already=already/earlyinf
           already2=already2/earlyinf
           faminfect=float(infectf)/earlyinf
           famnoinft=float(noinfectf)/earlyinf
      endif
      earlyR=0
      if(attac.gt.0) earlyR=earlyinf/(0.05*Nsize)
      if(suminf2.gt.0) faminfect2=float(infectf2)/suminf2
      write(72,185) faminfect,faminfect2,famnoinft,
     q               earlyR,earlyR*faminfect,earlyR*faminfect2,
     q               already*earlyR,already2*earlyR
 185  format('# Family infec/tot inf=',F6.3,'  When fam size=2:',F6.3,
     q       ' outside family=',F6.3,
     q       ' Re=',F6.2,', Prop infect in family=',2F5.2,
     q       ' fail:',2F5.2)
      write(72,178) 
      write(72,176) 
      write(72,173) int(sumc),(statgroup(i)/sumc,i=1,7)
      write(88,173) int(sumc),(statgroup(i)/sumc,i=1,7)
      write(72,183) int(sumx),(datamat(i,0,0)/sumx,i=1,6)
      write(71,173) 
      write(72,*)
 173  format('#',13x,I10,F9.3,6F10.3)
 183  format('# data to fit:',I10,F9.3,6F10.3)
 174  format('#',13x,F10.1,F9.4,6F10.5)
 179  format('#',17x,' par:    family  grandpar   crosinf',
     q       '      work    school:   others     w+s')
 176  format('#',17x,' stat.   family  grandpar   crosinf',
     q       '      work    school:   others     w+s')
 177  format('# activity in each of the sub-groups')
 178  format('# total fraction of activity in each of the sub-groups')
      do i=1,9
      write(71,171) i,statage(i),statage(i)/sum,polymod(i),attack(i)
      do j=1,9
      write(72,172) i,j,statcross(i,j)/sumc,
     q              statsub(1,i,j)/sumc,
     q              statsub(2,i,j)/sumc,
     q              statsub(3,i,j)/sumc,
     q              statsub(4,i,j)/sumc,
     q              statsub(5,i,j)/sumc
      enddo
      write(72,*)
      enddo
c
      write(73,*) 1,statgroup(1),faminf
      write(73,*) 2,statgroup(2),grandpar
      write(74,*) 3,statgroup(3),crosinf
      write(75,*) 4,statgroup(4),workinf
      write(76,*) 5,statgroup(5),unspec
      sumi0=0
      sumc0=0
      do i=0,Nsuperlim
         sumc0=sumc0+histacc(i)
         sumi0=sumi0+iacc(i)
      enddo
      sumi=0
      write(77,*) '#Nsuperlim=',Nsuperlim
      n80=0
      ssum=0
      sumi=0
      do i=0,Nsuperlim
      sumi=sumi+iacc(i)
      ssum=ssum+nncc(i)
      if(sumi.gt.0.2*sumi0.and.n80.eq.0) then
         n80=1
         ni=i
         ss=ssum
         ssi=1-sumi/sumi0
      endif
      enddo
c
      write(51,551) ni,1-ss/ssum,ssi,ssum,sumi,nncc(0) 
      write(72,551) ni,1-ss/ssum,ssi 
 551  format('# ni,s=',I6,2F9.4,2F12.1,I7)
      ymean=ymean/yaaa
      yvar=yvar/yaaa-ymean**2
      yvc=0
      if(yvar.gt.0) yvc=sqrt(yvar)/ymean
      write(177,778) ymean,yvar,yvc,yvc**2,(1./ymean+1./(1.-superpow))
 778  format('# mean R=',F7.2,' Var=',F10.2,
     q       ' cv,cv**2=',2F10.2,' expected=',2F10.3)
c
      sum=0
      sumn=0
      suma=0
      sumc=0
      yaaa=0
      ymean=0
      yvar=0
      do i=0,min(100,Nsuperlim)
      sum=sum+histsuper(i)
      sumc=sumc+histacc(i)
      sumn=sumn+nncc(i)
      suma=suma+float(iacc(i))/sumi
      yaaa=yaaa+nncc(i)
      ymean=ymean+i*nncc(i)
      yvar=yvar+i**2*nncc(i)
      anncc=nncc(i)/ssum
      if(i.gt.10) anncc=(nncc(i-1)+nncc(i)+nncc(i+1))/(3*ssum)
      write(177,777) i,histsuper(i),sum,
     q            anncc,float(iacc(i))/sumi,
     q            sumn/ssum,suma,histacc(i)
      enddo
      ymean=ymean/yaaa
      yvar=yvar/yaaa-ymean**2
      yvc=0
      if(yvar.gt.0.and.ymean.gt.0) yvc=sqrt(yvar)/ymean
      write(177,778) ymean,yvar,yvc,yvc**2,Rnow0,Vnow0
 777  format(' ',I6,6F13.8,2F12.4)
c
c
 187  format(' ',I5,5F12.6)
c
 100  format(' ',F8.1,3F8.4)
 120  format('# incubation=',F6.1,' serial time inf-inf=',F6.1,
     q       ' serial ons-ons',F6.1)
 130  format(' ',I2,I3,F7.3,2F7.4,I4,3F7.1,F9.3,2I4,3F8.4)
 170  format('#  samplet unspec=',F12.1,' sampled spec=',F12.1)
 171  format(' ',I6,F10.1,3F10.3)
 172  format(' ',2I6,6F10.5)
c
      stop
      end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0

      subroutine social(idum)
      common/sizes1/ worksize,famsize,frsize,crossize,schoolsize
      common/weights/ workinf,school,faminf,frinf,crosinf,ainf(9,9)
      common/weights1/ grandpar
      common/param1/ Nsize,konec,kmax
      common/society/ iage(1500000),iwork(200000,0:50),
     q                ifam(0:1000000,0:12),
     q     ifriend(0:500000,0:20)
      common/society1/ jwork(1500000),
     q            jfam(1500000),jfriend(1500000),jcross(1500000)
      common/society2/ Ncomp,Nschool
      common/ages/ jgroup(9,0:200000)
      common/links/ Kon(1100000),K(1100000,100)
      common/strength/ SK(1100000,100)
      common/links1/ Kiden(1100000,100)
      common/persons/ eta(1500000),agescal(10),actpers(1500000) 
      common/persons1/ actpmax,etamax,agemaxx,actpave,etaave
      common/persons2/ beta0,dt,actu,actp,unspec
      common/persons3/ agescal0(10)
      common/mossong/ polymod(9),polymodg(5)
      common/teachergroup/ iteacher(200000)
      dimension numage(200000)
      dimension iparents(1000000)
      dimension histfam(0:100),histchild(0:100)
      dimension histw(0:100),hists(0:100),histf(0:100),histc(0:100)
      dimension act(0:10,0:10)
      dimension histcoupl(0:100000),histpers(0:100000)
      iplot=0
      Nsize9=Nsize/9
      Nchildren=2*Nsize9
      Nactive=5*Nsize9
c
      Nfam4=Nsize/20
      Nfam2=Nsize/10
      Npair=Nsize/10
      Nsingle=Nsize-4*Nfam4-2*Nfam2-2*Npair
      ave=famsize+1     !parameter to set number of families
      avew=worksize+1    !Average company size (schools set to double)
      aves=schoolsize+1
c
c
      Numfam=Nsize/ave
c     jgroup(ig,x)= individ number x in age group ig, NOT USED YET
      do ig=1,9
          jgroup(ig,0)=0
      enddo
      do i=1,Nsize
         jfam(i)=0
         ig=1+(9*i)/Nsize
         if(ig.gt.9) ig=9
         iage(i)=ig
         jgroup(ig,0)=jgroup(ig,0)+1
         in=jgroup(ig,0)
         jgroup(ig,in+1)=i
      enddo
      do inum=1,Numfam
         ifam(inum,0)=0
         iparents(inum)=0
      enddo
c
      inum=0
      do i=1,Nsize
         if(jfam(i).eq.0) then
 10          inum=ran2(idum)*Numfam+1
c........    children only in 40% of families:
             if(iage(i).le.2) inum=0.4*ran2(idum)*Numfam+1
c........    no families with more than two parents
             if(iparents(inum).gt.2) goto 10
             if(ifam(inum,0).ge.11) goto 10
c........    a member of inum family:
             if(ifam(inum,0).gt.0.and.iage(i).gt.2) then
                 member=ifam(inum,1)
                 if(iage(member).lt.iage(i)-1.
     q               or.iage(member).gt.iage(i)+1) then
                     goto 10
                 endif
             endif
c........
             jfam(i)=inum                 !pointer from i to family number
             ifam(inum,0)=ifam(inum,0)+1 
             kkk=ifam(inum,0) 
             ifam(inum,kkk)=i             !pointer from family number to i
             if(iage(i).gt.2) then
                 iparents(inum)=iparents(inum)+1
                 ippp=iparents(inum)
             else if(iage(i).le.2) then
                 if(iparents(inum).le.2) then
                     ipar=i+(3*Nsize)/9
c
                     jfam(ipar)=inum                 !pointer from i to family number
                     ifam(inum,0)=ifam(inum,0)+1 
                     kkk=ifam(inum,0) 
                     ifam(inum,kkk)=ipar            !pointer from family number to i
                     iparents(inum)=iparents(inum)+1
                 endif
             endif
         endif
c........
      enddo
c
c.... Statistics of families (denmark average is 2.15):
c    
      numreal=0
      nnchild=0
      ak=0
      ac=0
      do ik=0,8
         histfam(ik)=0
         histchild(ik)=0
      enddo
      do inum=1,Numfam
         kkk=ifam(inum,0)
         if(kkk.gt.0) then
             numreal=numreal+1
             histfam(kkk)=histfam(kkk)+1
             ak=ak+kkk
             nchild=0
             do ik=1,kkk
                  kx=ifam(inum,ik)
                  kk=iage(kx)
                  if(kk.le.2) then
                      nchild=nchild+1
                  endif
             enddo
             nnchild=nnchild+nchild
             if(nchild.gt.0) then
             xc=xc+1
             ac=ac+nchild
             endif
             histchild(nchild)=histchild(nchild)+1
         endif
      enddo
      ak=ak/numreal
      write(14,91) Numfam,numreal,ak,nnchild,ac/xc
      do i=0,6
         write(14,*) i,histfam(i)/numreal,
     q               histchild(i)/numreal,ak/numreal,ac/xc
      enddo
c
c.....construction of interaction network with strength
c
      do i=1,Nsize
          Kon(i)=0
          do j=1,100
              Kiden(i,j)=0
              K(i,j)=0
              SK(i,j)=0
          enddo
      enddo
c
c.....Work part of life (for all below 70, and tripple size for children):
c
      Ncomp=Nactive/avew
      Nschool=Nchildren/aves
      Nschool2=Nschool/2
      Nteacher=Nschool*2
      write(81,*) Ncomp,Ncomp,Nsize,Nsize/9,avew
      do inum=1,Ncomp+Nschool
          iteacher(inum)=0
          do kk=0,50
          iwork(inum,kk)=0
          enddo
      enddo
      do i=1,Nsize
          jwork(i)=0
      enddo
      do i=1,Nsize
          ig=iage(i)
          if(ig.le.7) then
 20           if(ig.eq.1) then
                  inum=ran2(idum)*Nschool2+1
              else if(ig.eq.2) then
                  inum=ran2(idum)*Nschool2+Nschool2+1
              else if(ig.ge.3) then
                  inum=ran2(idum)*Ncomp+Nschool+1
              endif
              if(iwork(inum,0).ge.50) goto 20
              jwork(i)=inum
              iwork(inum,0)=iwork(inum,0)+1
              kk=iwork(inum,0)
              iwork(inum,kk)=i
          endif
      enddo
c
      do ik=0,30
         histw(ik)=0
         hists(ik)=0
      enddo
      aw=0
      do inum=1,Ncomp+Nschool
         kk=iwork(inum,0) 
         if(kk.gt.0) then
             k1=iwork(inum,1)
             ik=iwork(inum,kk)    !take one in work and use this as marker
             ig=iage(ik) 
             if(ig.gt.2) then
                 xw=xw+1
                 aw=aw+kk
                 histw(kk)=histw(kk)+1
             else
                 xs=xs+1
                 as=as+kk
                 hists(kk)=hists(kk)+1
             endif
         endif
      enddo
c
      do inumc=1,Nschool
 45       nsh=Nschool+1+ran2(idum)*Nteacher   !choose a teacher group
          iteacher(nsh)=iteacher(nsh)+1
          if(iteacher(nsh).ge.4) goto 45
          nth=iwork(nsh,0)                    !size of teacher group
          if(nth.le.0) goto 45
          do nxx=1,2                          !each class gets 2 teachers
 44       kt=nth*ran2(idum)+1                 !teacher to connect to
          itea=iwork(nsh,kt)                  !teacher indentifyer             
          if(itea.le.0) then
             goto 44
          endif
          nclass=iwork(inumc,0)               !class
          do i1=1,nclass                      !som over all children in group
               i=iwork(inumc,i1)              !child with identifyer
               if(Kon(itea).le.60) then
               Kon(i)=Kon(i)+1
               kk=Kon(i)
               K(i,kk)=itea
               Kiden(i,kk)=5
               ig=iage(i)
               factor=eta(i)*agescal(ig)
               SK(i,kk)=SK(i,kk)+school*factor
               Kon(itea)=Kon(itea)+1
               kk=Kon(itea)
               K(itea,kk)=i
               Kiden(i,kk)=5
               ig=iage(itea)
               factor=eta(itea)*agescal(ig)
               SK(itea,kk)=SK(itea,kk)+school*factor
               endif
          enddo
          enddo
      enddo
c
c
      aw=aw/xw
      as=as/xs
      write(15,92) Ncomp,xw,aw,Nschool,xs,as
      do i=0,30
         write(15,*) i,histw(i)/xw,hists(i)/xs
      enddo
c
c.....family networks:
c
      do i=1,Nsize
         jf=jfam(i)
         if(jf.gt.0) then
         ns=ifam(jf,0)
         if(ns.gt.1) then
         do is=1,ns
             in=ifam(jf,is)
             if(in.ne.i) then
                 Kon(i)=Kon(i)+1
                 kk=Kon(i)
                 if(kk.le.100) then
                     K(i,kk)=in
                     ig=iage(i)
                     factor=eta(i)*agescal(ig)
                     SK(i,kk)=SK(i,kk)+faminf*factor
                     Kiden(i,kk)=1
                 endif
             endif
         enddo
         endif
         endif
      enddo
      if(iplot.eq.0.or.iplot.eq.1) then
          write(18,*) '#grandpar+family'
          do i=1,25
          write(18,95) i,Kon(i),(K(i,j),j=1,10)
          enddo
      endif
c
c.....work networks:
      do i=1,Nsize
         jf=jwork(i)
         if(jf.gt.0) then
         ns=iwork(jf,0)
         if(ns.gt.0) then
         ig=iage(i)
         do is=1,ns
             if(iwork(jf,is).ne.i) then
                 Kon(i)=Kon(i)+1
                 kk=Kon(i)
                 if(kk.le.100) then
                 K(i,kk)=iwork(jf,is)
                 jn=k(i,kk)
                 jg=iage(jn)
                 facage=agescal(ig)
                 factor=eta(i)*agescal(ig)
                 if(ig.ge.3) SK(i,kk)=SK(i,kk)+workinf*factor
                 if(ig.ge.3) Kiden(i,kk)=4
                 if(ig.le.2) SK(i,kk)=SK(i,kk)+school*factor
                 if(ig.le.2) Kiden(i,kk)=5
                 endif
             endif
         enddo
         endif
         endif
      enddo
c
c.....projections:
c
      do i=1,9
      do j=1,9
          act(i,j)=0
      enddo
      enddo
      do i=0,100000
          histcoupl(i)=0
          histpers(i)=0
      enddo
      etamax=0
      agemaxx=0
      etaave=0
      do i=1,Nsize
          ia=iage(i)
          if(eta(i)*agescal(ia).gt.etamax) then
              etamax=eta(i)*agescal(ia)
          endif
          etaave=etaave+eta(i)*agescal(ia)
          if(agescal(ia).gt.agemaxx) then
              agemaxx=agescal(ia)
          endif
          ssum=0
          if(Kon(i).gt.0.and.Kon(i).le.100) then
          do ik=1,Kon(i)
             j=K(i,ik)
             if(j.gt.0.and.i.ne.j) then
                 jg=iage(j)
                 ssum=ssum+SK(i,ik)
                 act(ia,jg)=act(ia,jg)+SK(i,ik)
                 icoupl=10*SK(i,ik)
                 if(icoupl.lt.100000) then
                     histcoupl(icoupl)=histcoupl(icoupl)+1
                 endif
             endif
          enddo
          endif
          issum=10*ssum
          actpers(i)=ssum
          actpave=actpave+ssum
          if(issum.lt.100000) histpers(issum)=histpers(issum)+1
c         write(51,151) i,iage(i),eta(i),actpers(i),Kon(i),
c    q                  (SK(i,ki),ki=1,10)
      enddo
 151  format(' ',I5,I2,2F7.3,I4,10F8.4)
      etaave=etaave/Nsize
      actpave=actpave/Nsize
      actpmax=0
      do i=1,Nsize
         kk=Kon(i)
         do j=1,kk
            if(j.le.100) then
                SK(i,j)=SK(i,j)/actpers(i)
            else
                write(333,*) i,Kon(i),j,iage(i),iage(j)
                Kon(i)=100
            endif
         enddo
         if(actpers(i).gt.actpmax) actpmax=actpers(i)
      enddo
      call ageprofile(Nsize)
      do i=1,Nsize
          actpers(i)=actpers(i)/actpmax
      enddo
      do i=1,9
          write(19,96) (act(i,j)/Nsize,j=1,9)
          sum=0
          do j=1,9
          sum=sum+act(i,j)/Nsize
          write(29,*) i,j,act(i,j)/Nsize
          enddo
          write(28,*) i,sum,polymod(i)/30
          write(29,*)
      enddo
c
      do i=0,9999
         if(histcoupl(i).gt.0.or.histpers(i).gt.0) then
         write(37,*) float(i)/100.,histcoupl(i),histpers(i)
         endif
      enddo
      write(41,169) actpmax,etamax,agemaxx
      write(71,169) actpmax,etamax,agemaxx
c     
 81   format(' ',I6,':',20I6)
 91   format('# Numfam=',2I8,F8.2,I8,F8.2)
 92   format('#',I6,2F10.2,I6,2F10.2)
 93   format('#',I6,2F10.2)
 95   format(' ',2I4,20I5)
 96   format(' ',9F10.3)
 169  format('# maximal activities=',3F10.3)
      return
      end
      
      subroutine ageprofile(Nsize)
      common/mossong/ polymod(9),polymodg(5)
c     Data on age from Polymod:
      polymod(1)=12.5
      polymod(2)=18.
      polymod(3)=13.5
      polymod(4)=14.1
      polymod(5)=14.
      polymod(6)=12.3
      polymod(7)=9.2
      polymod(8)=7
      polymod(9)=7
      sum=0
      do i=1,9
          sum=sum+polymod(i)
      enddo
      do i=1,9
          polymod(i)=polymod(i)/sum
          write(28,*) i,sum,polymod(i)/30
      enddo
      return
      end

      subroutine initialize(idum)
      common/sizes1/ worksize,famsize,frsize,crossize,schoolsize
      common/weights/ workinf,school,faminf,frinf,crosinf,ainf(9,9)
      common/weights1/ grandpar
      common/param1/ Nsize,konec,kmax
      common/society/ iage(1500000),iwork(200000,0:50),
     q     ifam(0:1000000,0:12),
     q     ifriend(0:500000,0:20)
      common/society1/ jwork(1500000),
     q            jfam(1500000),jfriend(1500000),jcross(1500000)
      common/society2/ Ncomp,Nschool
      common/links/ Kon(1100000),K(1100000,100)
      common/links1/ Kiden(1100000,100)
      common/strength/ SK(1100000,100)
      common/ages/ jgroup(9,0:200000)
      common/persons/ eta(1500000),agescal(10),actpers(1500000) 
      common/persons1/ actpmax,etamax,agemaxx,actpave,etaave
      common/persons2/ beta0,dt,actu,actp,unspec
      common/persons3/ agescal0(10)
      common/mossong/ polymod(9),polymodg(5)
      common/hospitals/ hospital(10),death(10)
      common/hospitals1/ agecount(10),bedsicu(10)
      common/super1/ superfreq,superpow,cut
      common/super2/ Nsupermax,Nsuperlim,Nsuper(1500000)
      common/super4/ ssuper(1500000)
      common/bbc/ datamat(7,0:10,0:10)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c      
c.....interaction strength per link (divided with number in corresp. network)
      workinf=workinf/worksize
      school=school/schoolsize
      faminf=faminf/famsize  
      grandpar=grandpar/(2*famsize+2)
      crosinf=crosinf/crossize 
      netw=worksize+famsize+2*famsize+2+frsize+crossize
c     write(71,68) workinf,faminf,grandpar,crosinf,netw
c     write(51,68) workinf,faminf,grandpar,crosinf,netw
 68   format('# couplings=',4F9.4,I6)
      do i=1,Nsize
         Nsuper(i)=1
         ssuper(i)=1
         if(ran2(idum).lt.superfreq) then
             Nsuper(i)=Nsupermax
             ssuper(i)=Nsupermax
         endif
      enddo
c
c.....Social activity as function of age: (adjusted to Mossong et al) 
c
      do i=1,10
         agescal(i)=1
      enddo
c
c.....agescal must be below 2:
      agescal(1)=0.88  !0.5
      agescal(2)=1.18
      agescal(3)=1.55   !1.6
      agescal(4)=1.1 
      agescal(5)=1.05 
      agescal(6)=1.4   !not so much fromm families, but same overall activity
      agescal(7)=0.98
      agescal(8)=1.2    !1.04
      agescal(9)=1.2   !1.14
c.....so that all agescal number less than 1 
c.....(can draw random number and compare):
      do i=1,10
         agescal(i)=agescal(i)/2
         agescal0(i)=agescal(i)
      enddo
c     option for heterogeneous social activity:
      do i=1,Nsize
         eta(i)=-beta0*log(ran2(idum))
         eta(i)=beta0
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Other option for coupling, (old)
c     only need subroutine social now
c
      do i=1,Nsize
          Kon(i)=0
          do j=1,100
             K(i,j)=0
             Kiden(i,j)=0
             SK(i,j)=0
          enddo
      enddo
c
      power=0.5
      kmax=100
      call social(idum)
c
c     faminf0,grandpar0,crosinf0,workinf0,unspec0,school0,crosinf0
      home=2*(22+62+40)+31+31+42
      work=2*(13+25+7)+2+35+6
      sch=54+4+1+2*(9+10+1)
c
      other=12+30+51+2*(16+43+51)
      sss=home+work+sch+other
      write(72,272) home/sss,work/sss,sch/sss,other/sss
 272  format('# home,work,sch,other=',4F8.3) 
      datamat(1,0,0)=home
c other in BBC:  2*(16+43+51)+12+30+51-->grandpar+friends+unspec
      datamat(2,0,0)=(16+43+51)/2           !grandpar (from other)
      datamat(3,0,0)=(12+30+51)             !crosinf (from Other)
      datamat(4,0,0)=work                   !work
      datamat(5,0,0)=sch                    !school
      datamat(6,0,0)=16+43+51+(16+43+51)/2  !unspec, from other
      do i=1,9
      do j=1,9
          datamat(1,i,j)=0.
      enddo
      enddo
      datamat(1,1,1)=31./4
      datamat(1,1,2)=31./4
      datamat(1,2,1)=31./4.
      datamat(1,2,2)=31./4.
c
      do i=1,2
      do j=3,6
c.........from i to j: (from top row to column)
          datamat(1,i,j)=62./8.
          datamat(1,j,i)=62./8.
      enddo
      enddo
      do i=3,6
          datamat(1,i,i)=31./4.   !smae age mariages.
      enddo
c
      do i=7,9
          do j=1,2
c.........from i to j: (from top row to column)
              datamat(1,i,j)=40./4.
          enddo
          do j=3,6
c.........from i to j: (from top row to column)
          datamat(1,i,j)=22./8.
          enddo
          datamat(1,i,i)=42./2.
      enddo
      do j=7,9
          do i=1,2
              datamat(1,i,j)=40./4.
          enddo
          do i=3,6
          datamat(1,i,j)=22./8.
          enddo
      enddo
c
c.....correction factor from equal many in each group, to Daish age profile:
      agecount(1)=626506
      agecount(2)=684412
      agecount(3)=765439
      agecount(4)=673969
      agecount(5)=781712
      agecount(6)=781703
      agecount(7)=672152
      agecount(8)=513155
      agecount(9)=249721
      sum=0
      do i=1,9
         sum=sum+agecount(i)
      enddo
      do i=1,9
         agecount(i)=agecount(i)/sum
         write(69,*) i,agecount(i),sum
      enddo
c
c.....propability to end in hspital, given age group
      hospital(1)=0
      hospital(2)=0.00013
      hospital(3)=0.0037
      hospital(4)=0.014
      hospital(5)=0.027
      hospital(6)=0.027
      hospital(7)=0.039
      hospital(8)=0.055
      hospital(9)=0.055
c.....propability to end at ICU given in hospital
      bedsicu(1)=0.05
      bedsicu(2)=0.05
      bedsicu(3)=0.05
      bedsicu(4)=0.05
      bedsicu(5)=0.063
      bedsicu(6)=0.122
      bedsicu(7)=0.274
      bedsicu(8)=0.432
      bedsicu(9)=0.709
      do i=1,9
         bedsicu(i)=bedsicu(i)*hospital(i)
      enddo
c.....spanish death rates converted to overall IFR fatality of 0.003
c.....danish serum tests suggest instead this should be        0.00232
      death(1)=0.00029
      death(2)=0.00029
      death(3)=0.00020
      death(4)=0.00029
      death(5)=0.00059
      death(6)=0.00137
      death(7)=0.00450
      death(8)=0.01331
      death(9)=0.02055
c
      return
      end

      subroutine sspower(idum)
      common/param1/ Nsize,konec,kmax
      common/super1/ superfreq,superpow,cut
      common/super2/ Nsupermax,Nsuperlim,Nsuper(1500000)
      common/super3/ histsuper(0:10000)
      common/super4/ ssuper(1500000)
      dimension pss(10000),hist(0:10000),gini(0:10000) 
      dimension spss(10000)
      dimension accinf(0:10000)
      dimension  gini0(0:1000)
      write(333,*) superpow,Nsuperlim
      disp=1-superpow
      dx=0.1
      pss(1)=1.
      pp=dx/2
      xx=1
      aa=1
c calculate right infection in first bin
      if(superpow.lt.1) then
          x=dx
          if(superpow.gt.0) then
             pss(1)=x**(1-superpow)*exp(-x/cut)/(1-superpow)
             spss(1)=x**(2-superpow)*exp(-x/cut)/(2-superpow)
             pp=spss(1)/pss(1)
          else
             pss(1)=exp(-x/cut)
             spss(1)=pp*exp(-x/cut)
          endif
          xx=pss(1)
          aa=(x/2)*pss(1)
      endif
      do i=2,Nsuperlim
         x=dx*float(i)-dx/2
         dummy=exp(-x/cut)*dx
         dummys=x*exp(-x/cut)*dx
         if(superpow.gt.1) then
             x=i
             dummy=x**(-superpow)
         else if(superpow.lt.1.and.superpow.gt.0) then
             dummy=x**(-superpow)*exp(-x/cut)*dx
             dummys=x*dummy
         endif
         pss(i)=pss(i-1)+dummy
         spss(i)=spss(i-1)+dummys
         xx=xx+dummy
         aa=aa+x*dummy
      enddo
      do i=0,Nsuperlim
           hist(i)=0
      enddo
      aa=aa/xx
      do i=1,Nsuperlim
         pss(i)=pss(i)/pss(Nsuperlim)
         spss(i)=spss(i)/spss(Nsuperlim)
      enddo
      do i=1,Nsize
         xx=ran2(idum)
         do ic=1,Nsuperlim
            if(pss(ic).gt.xx) then
                isuper=ic
                goto 111
            endif
         enddo
 111     Nsuper(i)=dx*isuper
         if(Nsuper(i).le.1) Nsuper(i)=1
         ssuper(i)=dx*(isuper-0.5)
         if(isuper.eq.1) ssuper(i)=pp
         hist(isuper)=hist(isuper)+1
      enddo
      accinf(0)=0
      ave=0
      var=0
      do i=1,Nsuperlim
         hist(i)=hist(i)/Nsize
         accinf(i)=accinf(i-1)+dx*(i-0.5)*hist(i)
         ave=ave+(dx*i-dx/2)*hist(i) 
         var=var+(dx*i-dx/2)**2*hist(i) 
      enddo
      var=var-ave**2
      spread=sqrt(var)
      cv=spread/ave
      do i=1,Nsuperlim
         accinf(i)=accinf(i)/accinf(Nsuperlim)
      enddo
      i80=0
      sum=0
      do i=1,Nsuperlim
         sum=sum+hist(i)
         histsuper(i)=hist(i)
         if(accinf(i).gt.0.2.and.i80.eq.0) then
            i80=1
            n80=i
            pop80=1-sum
            act80=1-accinf(i)
         endif
         write(77,*) i*dx,pss(i),spss(i),
     q               histsuper(i),sum,accinf(i)
      enddo
      analytical=(1-superpow)*cut
      write(77,877) Nsuperlim,cut,disp,dx,n80,pop80,
     q        act80,ave,spread,cv
      write(51,877) Nsuperlim,cut,disp,dx,n80,pop80,
     q        act80,analytical,aa,ave,spread,cv
      write(72,877) Nsuperlim,cut,disp,dx,n80,pop80,
     q        act80,analytical,aa,ave,spread,cv
 877  format('# Nlim,cut,k=',I6,F7.1,F7.3,F7.2,
     q       ' n80,p,a=',I4,2F9.3,' anal=',2F7.2,
     q       ' ave,spread,cv=',3F9.2)
      return
      end

      subroutine data
      common/hospitalization/ lockdown1,inhosp(-300:1000)
      common/sweeden/ lockdowns1,insw(-300:1000),ideaths(-300:1000)
      dimension smooth(-300:1000)
      lockdown1=18
      lockdowns1=33   !swedish lock down time
      lockdowns1=31   !swedish lock down time
      lockdowns1=29   !swedish lock down time
      do i=-300,1000
          ideaths(i)=-401
          inhosp(i)=-401
          insw(i)=-401
      enddo
c
      inhosp(-3)=0
      inhosp(-2)=0
      inhosp(-1)=0
      inhosp(0)=0
      inhosp(1)=0    !25/2 2020  
      inhosp(2)=0      
      inhosp(3)=0      
      inhosp(4)=0      
      inhosp(5)=0      
      inhosp(6)=1
      inhosp(7)=0
      inhosp(8)=1
      inhosp(9)=0  
      inhosp(10)=1   
      inhosp(11)=4
      inhosp(12)=0
      inhosp(13)=3 
      inhosp(14)=9
      inhosp(15)=15
      inhosp(16)=10
      inhosp(17)=18
      inhosp(18)=26
      inhosp(19)=24
      inhosp(20)=33
      inhosp(21)=45
      inhosp(22)=47
      inhosp(23)=47
      inhosp(24)=52
      inhosp(25)=64
      inhosp(26)=24
      inhosp(27)=40
      inhosp(28)=86
      inhosp(29)=85
      inhosp(30)=94
      inhosp(31)=92
      inhosp(32)=90
      inhosp(33)=91
      inhosp(34)=51
      inhosp(35)=84
      inhosp(36)=57
      inhosp(37)=61
      inhosp(38)=73
      inhosp(39)=72
      inhosp(40)=45
      inhosp(41)=34
      inhosp(42)=52
      inhosp(43)=56
      inhosp(44)=28
      inhosp(45)=34
      inhosp(46)=34
      inhosp(47)=31
      inhosp(48)=27
      inhosp(49)=30
      inhosp(50)=46
      inhosp(51)=30
      inhosp(52)=33
      inhosp(53)=27
      inhosp(54)=27
      inhosp(55)=34
      inhosp(56)=33
      inhosp(57)=38
      inhosp(58)=28
      inhosp(59)=33
      inhosp(60)=23
      inhosp(61)=27
      inhosp(62)=26
      inhosp(63)=21
      inhosp(64)=15
      inhosp(65)=15
      inhosp(66)=16
      inhosp(67)=16
      inhosp(68)=24
      inhosp(69)=16
      inhosp(70)=15
      inhosp(71)=16
      inhosp(72)=6
      inhosp(73)=15
      inhosp(74)=9
      inhosp(75)=10
      inhosp(76)=10
      inhosp(77)=8
      inhosp(78)=7
      inhosp(79)=6
      inhosp(80)=13
      inhosp(81)=17
      inhosp(82)=5
      inhosp(83)=4
      inhosp(84)=14
      inhosp(85)=7
      inhosp(86)=7
      inhosp(87)=6
      inhosp(88)=13
      inhosp(89)=9
      inhosp(90)=2
      inhosp(91)=6
      inhosp(92)=11
      inhosp(93)=7
      inhosp(94)=5
      inhosp(95)=9
      inhosp(96)=7
      inhosp(97)=4
      inhosp(98)=5
      inhosp(99)=6
      inhosp(100)=3
      inhosp(101)=5
      inhosp(102)=2
      inhosp(103)=4
      inhosp(104)=3
c
c
      insw(-5)=0
      insw(-4)=0
      insw(-3)=0
      insw(-2)=0
      insw(-1)=0
      insw(0)=0
      insw(1)=0   !25/2
      insw(2)=1
      insw(3)=1
      insw(4)=8
      insw(5)=3
      insw(6)=0
      insw(7)=5
      insw(8)=13
      insw(9)=30
      insw(10)=25
      insw(11)=59
      insw(12)=33
      insw(13)=46
      insw(14)=101
      insw(15)=98
      insw(16)=196
      insw(17)=151
      insw(18)=152
      insw(19)=71
      insw(20)=69   !15/3
      insw(21)=83
      insw(22)=119
      insw(23)=145
      insw(24)=143
      insw(25)=180
      insw(26)=134
      insw(27)=117
      insw(28)=117   !fail
      insw(29)=182
      insw(30)=230
      insw(31)=314
      insw(32)=286
      insw(33)=366
      insw(34)=300
      insw(35)=281
      insw(36)=416
      insw(37)=475  !31/3
      insw(38)=486
      insw(39)=554
      insw(40)=601
      insw(41)=357
      insw(42)=340
      insw(43)=389
      insw(44)=738
      insw(45)=655
      insw(46)=645
      insw(47)=454
      insw(48)=395
      insw(49)=464
      insw(50)=437
      insw(51)=480
      insw(52)=604
      insw(53)=623
      insw(54)=688
      insw(55)=532
      insw(56)=388
      insw(57)=461
      insw(58)=706
      insw(59)=722
      insw(60)=758
      insw(61)=786
      insw(62)=473   !
      insw(63)=300
      insw(64)=576
      insw(65)=761
      insw(66)=830
      insw(67)=678
      insw(68)=532   !1 st may
      insw(69)=298
      insw(70)=262
      insw(70)=501
      insw(71)=662
      insw(72)=751
      insw(73)=821
      insw(74)=713
      insw(75)=509
      insw(76)=278
      insw(77)=486
      insw(78)=800
      insw(79)=721
      insw(80)=686
      insw(81)=711
      insw(82)=358
      insw(83)=259
      insw(84)=457
      insw(85)=689
      insw(86)=825
      insw(87)=614   !21 may 2020
      insw(88)=529   !
      insw(89)=403   !
      insw(90)=210   !
      insw(91)=505   !
      insw(92)=505   !
      insw(93)=757   !
      insw(94)=816   !
      insw(95)=771   !
      insw(96)=776   !
      insw(97)=432   !
      insw(98)=265   !1 juni
      insw(99)=644   !
      insw(100)=889  !
      insw(101)=1070 !
      insw(102)=1042 !
      insw(103)=1016 !
      insw(103)=500 !
      insw(103)=500 !
      insw(103)=1000!
c
c     deaths in whole Sweden:
      ideaths(1)=0    !
      ideaths(2)=0    !
      ideaths(3)=0    !
      ideaths(4)=0    !
      ideaths(5)=0    !
      ideaths(6)=0    !
      ideaths(7)=0    !
      ideaths(8)=0    !
      ideaths(9)=0    !
      ideaths(11)=0    !
      ideaths(12)=0    !
      ideaths(13)=0    !
      ideaths(14)=0    !
      ideaths(15)=0    !
      ideaths(16)=1    !11 march
      ideaths(17)=0    !12 march
      ideaths(18)=1    !13 march
      ideaths(19)=1    !14 march
      ideaths(20)=2    !15 march
      ideaths(21)=2    !16 march
      ideaths(22)=1    !17 march
      ideaths(23)=6    !18 march
      ideaths(24)=7    !19 march
      ideaths(25)=9    !20 march
      ideaths(26)=8    !21 march
      ideaths(27)=11   !22 march
      ideaths(28)=11   !23 march
      ideaths(29)=21   !24 march
      ideaths(30)=22   !25 march
      ideaths(31)=31   !26 march
      ideaths(32)=32   !27 march
      ideaths(33)=35   !28 march
      ideaths(34)=38   !29 march
      ideaths(35)=45   !30 march
      ideaths(36)=48   !31 march
      ideaths(37)=53   ! 1 april
      ideaths(38)=70   ! 2 april
      ideaths(39)=79   !
      ideaths(40)=70   !
      ideaths(41)=86   !
      ideaths(42)=90   !
      ideaths(43)=84   !
      ideaths(44)=115  ! 8 april
      ideaths(45)=86   !
      ideaths(46)=90   !
      ideaths(47)=103  !
      ideaths(48)=97   !
      ideaths(49)=85   !
      ideaths(50)=91   !14 april
      ideaths(51)=115  !
      ideaths(52)=111  !
      ideaths(53)=82   !
      ideaths(54)=86   !
      ideaths(55)=88   !
      ideaths(56)=85   !
      ideaths(57)=62   !
      ideaths(58)=77   !
      ideaths(59)=86   !23 april
      ideaths(60)=89   !
      ideaths(61)=73   !
      ideaths(62)=75   !
      ideaths(63)=73   !
      ideaths(64)=82   !
      ideaths(65)=83   !
      ideaths(66)=78   !
      ideaths(67)=78   !
      ideaths(68)=73   !
      ideaths(69)=75   !
      ideaths(70)=84   ! 4 may
      ideaths(71)=72   !
      ideaths(72)=73   !
      ideaths(73)=80   !
      ideaths(74)=60   !
      ideaths(75)=67   !
      ideaths(76)=74   !
      ideaths(77)=64   !
      ideaths(78)=60   !12 may
      ideaths(79)=51   !
      ideaths(80)=46   !14 may
      ideaths(81)=57   !
      ideaths(82)=48   !
      ideaths(83)=53   !
      ideaths(84)=61   !
      ideaths(85)=39   !
      ideaths(86)=53   !
      ideaths(87)=52   !
      ideaths(88)=54   !
      ideaths(89)=56   !
      ideaths(90)=42   !24 may
      ideaths(91)=41   !
      ideaths(92)=27   !
      ideaths(93)=39   !
      ideaths(94)=38   !
      ideaths(95)=39   !
      ideaths(96)=37   !30 may
      ideaths(97)=43   !
      ideaths(98)=39   ! 1 june
      ideaths(99)=29   !
      ideaths(100)=22  !
      ideaths(101)=40  !
      ideaths(102)=25  !
      ideaths(103)=13  !
      ideaths(104)=18  !
      ideaths(105)=21  !
      ideaths(106)=10  !
      ideaths(107)=12  ! 10 june
      ideaths(108)=8  !
      ideaths(109)=5  !
      ideaths(108)=-401
      ideaths(109)=-401
c     stockholm region instead:
      insw(-5)=0
      insw(-4)=0
      insw(-3)=0
      insw(-2)=0
      insw(-1)=0
      insw(0)=0
      insw(1)=0   !25/2
      insw(2)=1
      insw(3)=1
      insw(4)=2
      insw(5)=1
      insw(6)=0
      insw(7)=1
      insw(8)=10
      insw(9)=21
      insw(10)=22
      insw(11)=36
      insw(12)=21
      insw(13)=29
      insw(14)=64
      insw(15)=26
      insw(16)=32
      insw(17)=42
      insw(18)=31
      insw(19)=18
      insw(20)=17   !15/3
      insw(21)=34
      insw(22)=35
      insw(23)=58
      insw(24)=66
      insw(25)=84
      insw(26)=71
      insw(27)=59
      insw(28)=99
      insw(29)=105
      insw(30)=154
      insw(31)=132
      insw(32)=176
      insw(33)=147
      insw(34)=150
      insw(35)=172
      insw(36)=209  !31/3
      insw(37)=205
      insw(38)=216
      insw(39)=245
      insw(40)=129
      insw(41)=172
      insw(42)=131
      insw(43)=243
      insw(44)=271
      insw(45)=240
      insw(46)=148
      insw(47)=200
      insw(48)=182
      insw(49)=200
      insw(50)=179
      insw(51)=215
      insw(52)=221
      insw(53)=221
      insw(54)=180
      insw(55)=192
      insw(56)=211
      insw(57)=163
      insw(58)=288
      insw(59)=291
      insw(60)=233
      insw(61)=138   !
      insw(62)=110
      insw(63)=226
      insw(64)=259
      insw(65)=279
      insw(66)=257
      insw(67)=141   !1 st may
      insw(68)=80
      insw(69)=128
      insw(70)=173
      insw(71)=173
      insw(72)=212
      insw(73)=304
      insw(74)=235
      insw(75)=114
      insw(76)=78
      insw(77)=241
      insw(78)=280
      insw(79)=259
      insw(80)=177
      insw(81)=207
      insw(82)=64
      insw(83)=58
      insw(84)=177
      insw(85)=203
      insw(86)=191
      insw(87)=124   !21 may 2020
      insw(88)=146   !
      insw(89)=78   !
      insw(90)=59   !
      insw(91)=181   !
      insw(92)=213   !
      insw(93)=225   !
      insw(94)=282   !
      insw(95)=215   !
      insw(96)=73   !
      insw(97)=67   !
      insw(98)=160   !1 juni
      insw(99)=152   !
      insw(100)=250  !
      insw(101)=219 !
      insw(102)=235 !
      insw(103)=80 !false (27,114)
      insw(104)=80 !false
      insw(105)=230 !
      insw(106)=230 !
      do i=1,100
          sum=0
          index=i-2
          do jj=1,5
            j=index+jj 
            sum=sum+insw(j)
          enddo
          smooth(i)=sum/5 
      enddo
c
      do i=1,100
          write(777,*) i,insw(i),smooth(i)
          insw(i)=smooth(i)
      enddo
c
      return
      end


