svp <- function (x){    #function to convert saturation vapour pressure units
    svp = 6.11*10^((7.5*x)/(237.7+x))    #saturation vapour pressure (hpa)
    #svp_m = 0.750061683 * svp             #svp in  mmhg
    return (svp/10)  #hPa to kPa
}

s2nk<-function(N,exp,phase,interval){
  #find out which exp and phase
  case=paste(exp,phase,sep="_")
  shift=0
  start=switch(case,
               H_1 = as.POSIXlt("2014-09-02 14:35:00"),
               H_2 = as.POSIXlt("2014-09-02 16:31:00"),        
               K_1 = as.POSIXlt("2014-09-04 14:35:00"),
               K_2 = as.POSIXlt("2014-09-04 16:33:00"),           
               N_1 = as.POSIXlt("2014-09-16 14:40:20"),
               N_2 = as.POSIXlt("2014-09-16 16:25:00"),
               
  )
  return(
    subset(
      get(paste(N,"_",exp,sep=""))
      ,FORMATTED.DATE.TIME >= start+shift & FORMATTED.DATE.TIME <= start+shift + interval*60)
  )
}

sub <- function(data){

    # Initial temperatures
    tcl = tclold = 0     # To prevent while loop from bugging out
    tcr = 36.9
    tsk = 34
    ttsk = 33.7     # setpoint tsk
    ttcr = 36.8     # setpoint tcr
    ata = 1 # atmospheres?
  
    # Clothing related
    clo = 0.57
    chclo = 1/(0.155*clo)    # pants + t-shirt
    
    facl = 1.0+0.15*clo     # surface enlargement?

    # Vascular
    skbf = 6.3
    skbfl = 90 # L/m2/hr

    # Sweat
    regswl = 500
    csw = 170 # g/m2/hr
    cdil = 75 #litres/(m2/h/K)  #200 in 1986
    cstr = 0.5                  #0.1 in 1986
    sweat = 0

    # Metabolism and activity
    mets = 2.9
    rm = 58.2
    me = 0.2
    
    # Body metrics
    ht = 1.81
    wt = 70
    adu = 1.8
    alpha = 0.044 + 0.35/ (skbf-0.1386)

    esk = 7.3 #init
    bz = 0.1

    # Time
    time = 0.0      # init
    exp_time = 30.0 # mins
    interval = 20   # seconds
    steps = exp_time*60/interval

    # Variable vectors
    vmrt = c()
    vtsk = c()
    vdtsk = c()
    vtcr = c()
    vtcl = c()
    vr = c()
    vc = c()
    ve = c()
    vedif = c()
    versw = c()
    vhfsk = c()
    vhfcr = c()
    

    index = 0
    while (time < 1.0){
        if (time >= 0.5) {
          act = rm *mets #no mets
        }
        else {
          act = rm * mets
        }
        index=index+1
        print(paste(index,"---------------"))
        time = time + 20/1800                
        
        #dynamic weather
        ta = data$TP[index]
        tr = data$MRT[index]
        v = data$WS[index]
        rh = data$RH[index]/100
        pa = rh * svp(ta)
        #####
        # Radiative and convective
        ####
        
        #still air
        chc1 = 5.66*(mets-0.85)^0.39  # based on walking activity in still air
        chc2=8.6*(v)^0.5
        if (time <= 0.5){
          chc=max(chc2,chc1)
        }
        else{
          chc=max(chc2,3)
        
        }
        
        chr = 4*0.72*5.67*10^-8*((tcl+tr)/2+273.15)^3
        
        tcl = (chclo*tsk+facl*(chc*ta+chr*tr))/(chclo+facl*(chc+chr))
        
        while (abs(tcl-tclold) > 0.01){
            tclold=tcl
            chr = 4*0.72*5.67*10^-8*((tcl+tr)/2+273.15)^3
            tcl = (chclo*tsk+facl*(chc*ta+chr*tr))/(chclo+facl*(chc+chr))           
        }        
        
        # heat flow from clothing to environment
        rad = facl*chr*(tcl-tr)
        conv = facl*chc*(tcl-ta)
        dry = rad+conv   # R + C negative = gain
        #dry = facl*(chc*(tcl-ta)+chr*(tcl-tr))   # R + C negative = gain

        # dry and latent respiratory heat losses
        eres = 0.017251 * act * (5.8662-pa)      # evaporative loss from respiration
        #eres2=0.0023*rm*(44-rh*svp(ta))  #1971 formulation

        cres=0.0014 * act * (34.0-ta) # *ata*ff  # convective loss from respiration
        
        wk = act * me
        # heat flow from skin
        hfsk=(tcr-tsk)*(5.28+1.163*skbf)-dry-esk
        
        # heat flow from core
        hfcr=act-(tcr-tsk)*(5.28+1.163*skbf)-cres-eres-wk

        # thermal capacities
        tccr = 58.2*(1-alpha)*wt *20   # 58.2 = 0.97 * 60s!!!!
        tcsk = 58.2*alpha*wt *20
    
        # temperature change in 1 min
        dtsk = (hfsk*adu)/tcsk 
        dtcr = (hfcr*adu)/tccr 
        dtbm=alpha*dtsk+(1.-alpha)*dtcr
        
        tsk = tsk+dtsk
        tcr = tcr+dtcr

       
        if (tsk>ttsk) {
          warms = tsk - ttsk
          colds = 0
        }
        else {
          colds = ttsk - tsk
          warms = 0
        }

        if (tcr > ttcr) {
          warmc = tcr - ttcr
          coldc = 0
        }
        else{
          coldc = ttcr - tcr
          warmc = 0
        }

        ttbm = bz*ttsk+(1-bz)*ttcr
        tbm = alpha*tsk + (1-alpha)*tcr
        if (tbm > ttbm){
          warmb = tbm - ttbm
          coldb = 0
        }
        else{
          coldb = ttbm - tbm
          warmb = 0
        }

        # temperature regulation

        dilat = cdil*warmc
        stric = cstr*colds
        
        skbf = (6.3+dilat)/(1+stric)
        # limits
        if (skbf < 0.5) {skbf = 0.5}
        if (skbf > skbfl) {skbf = skbfl}

        # skin-core proportion changes with blood flow
        alpha = 0.0417737+ 0.7451832/(skbf+0.58517)
        #alpha = 0.1

#         cat('warms',warms,'\n')
        #sweating regulation
        regsw = csw*warmb*exp(warms/10.7)
        if (regsw > regswl) {
          regsw = regswl
        }
        
        ersw = 0.68 * regsw
#        ersw = 0.7 * regsw * 2^((tsk-ttsk)/3)   #1971
        #shivering
        actold = act
        act = act + 19.4 * colds*coldc
        
        mshiv = 19.4 * colds*coldc
       
        #evaporation
        lr=15.1512*(tsk+273.15)/273.15
        im = 1
        icl = 0.45
        rt = (1/im) * (1/(lr*facl*chc)+1/(lr*chclo*icl))        
        emax = (1/rt)*(svp(tsk)-pa)
        
        prsw = ersw / emax

        pdif = (1-prsw)*0.06
        edif = pdif*emax
        esk = ersw + edif
        sweat = sweat + (esk*adu/0.7)*(interval/3600) #g/timestep
        pwet = esk/emax

        #dripping sweat
        eveff = 1
        
        if (pwet >= eveff & emax >= 0){
            pwet = eveff
            prsw = (eveff-0.06)/.94
            ersw = prsw / emax
            pdif = (1-prsw)*0.06
            edif = pdif*emax
            esk = ersw + edif
        }

        if (emax < 0) {
            pdif = 0
            edif = 0
            esk = emax
            pwet = eveff
            prsw = eveff
            ersw= 0
        }
        edrip = (regsw*0.68-prsw*emax)/0.68
        if (edrip < 0) {edrip = 0}

        #vapour pressure at skin
        vpsk = pwet * svp(tsk)+(1-pwet)*pa
        # rh at skin
        rhsk = vpsk/svp(tsk)
#         cat("tsk:",tsk,"\n")
        vtsk <- c(vtsk,tsk)
        vdtsk <- c(vdtsk,dtsk)
        vtcr <- c(vtcr,tcr)
        vtcl <- c(vtcl,tcl)
        vr <- c(vr,rad)
        vc <- c(vc,conv)
        ve <- c(ve,esk)
        vmrt <- c(vmrt,tr)
        vhfsk <- c(vhfsk,hfsk)
        vhfcr <- c(vhfcr,hfcr)
        vedif <- c(vedif,edif)
        versw <- c(versw,ersw)

        }
        par(xaxs='i',yaxs='i',xpd=F,cex=1)
        old.par <- par(mfrow=c(3,2),mar=c(2,2,2,2))
        plot(vhfsk,type="l",main="hfsk and hfcr (red)",ylim=c(-100,100))
        lines(vhfcr,col="red")
        plot(vtsk,type="l",ylim=c(30,36), main= "tsk and tcl (red)")
        lines(vtcl,col="red")
        plot(vtcr,type="l",ylim=c(36.5,37.5))
        plot(vr,type="l", main="R")
        plot(vc,type="l", main="C")
        plot(ve,type="l", main = "E",ylim=c(0,120))
        lines(vedif,col="blue")
        lines(versw,col="red")
        
        cat('sweat',sweat,'\n')
        return(list(tsk=vtsk,tcr=vtcr))
    }

  data<-s2nk("K01","K",1,30)
  sub(data)


sub2 <- function(){
  
  ctc = chc + chr
  to = ta + erf / ctc
  cloe = clo - (facl - 1)/(0.155 * facl * ctc)
  fcle = 1/(1 + 0.155 * ctc * cloe)
  fpcl = 1/(1 + (0.155/icl) * chc * cloe )
  
  hsk = ctc * fcle * (tsk - to) + pwet * lr * chc * fpcl * (svp(tsk) - pa)
  
  et = tsk - hsk / (ctc * fcle)
  
  # approximate ET*
  while (err >= 0) {
    et = et + 0.1
    err = hsk - ctc * fcle * (tsk - et) - pwet * lr * chc * fpcl * (svp(tsk) - svp(et)/2)
  }
  
  # SET*
  chrs = chr
  chcs = 5.66 * (act/58.2-0.85)^0.39
  if (chcs < 3) {
    chcs = 3
  }
  
  rn = rm -wk
  clos = 1.3264/(rn/58.15+0.7383)-0.0953
  kclos = 0.25
  facls = 1+kclos + clos
  ctcs = chrs + chcs
  cloes = clos - (facls -1)/(0.155*facls*ctcs)
  fcles = 1/(1+0.155*ctcs*cloes)
  fpcls = 1/(1+(0.155/0.45)*chcs*cloes)
  
}
