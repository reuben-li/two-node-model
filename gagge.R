svp <- function (x){    #function to convert saturation vapour pressure units
    svp = 6.11*10^((7.5*x)/(237.7+x))    #saturation vapour pressure (hpa)
    return (svp/10)  #hPa to kPa
}

svp_m <- function (x){    #function to convert saturation vapour pressure units
  svp = 6.11*10^((7.5*x)/(237.7+x))    #saturation vapour pressure (hpa)
  svp_m = 0.750061683 * svp             #svp in  mmhg
  return (svp_m)  #hPa to kPa
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
               N_1 = as.POSIXlt("2014-09-16 14:42:20"),
               N_2 = as.POSIXlt("2014-09-16 16:26:00"),
               
  )
  return(
    subset(
      get(paste(N,"_",exp,sep=""))
      ,FORMATTED.DATE.TIME >= start+shift & FORMATTED.DATE.TIME <= start+shift + interval*60)
  )
}

s2n<-function(N,exp,phase,interval){
  #find out which exp and phase
  case=paste(exp,phase,sep="_")
  shift=0
  start=switch(case,
               H_1 = as.POSIXlt("2014-09-02 14:35:00"),
               H_2 = as.POSIXlt("2014-09-02 16:31:00"),        
               K_1 = as.POSIXlt("2014-09-04 14:35:00"),
               K_2 = as.POSIXlt("2014-09-04 16:33:00"),           
               N_1 = as.POSIXlt("2014-09-16 14:38:40"),
               N_2 = as.POSIXlt("2014-09-16 16:25:00"),
               
  )
  return(
    subset(
      get(paste(N,"_",exp,sep=""))
      ,TIME > start+shift & TIME < start + interval*60+shift)
  )
}

sub <- function(data,phys,exp,n,phase){

    # Initial temperatures
    tcl = tclold = 0     # To prevent while loop from bugging out
    tcr = 36.9
    tcr = phys$ttcr # use base oral temp as start temp # higher figure
#     tcr=s2n(paste("N0",n,sep=""),exp,phase,30)$Ch.1[1]
    tsk = 34
    tsk=s2n(paste("N0",n+4,sep=""),exp,phase,30)$mean[1]
#     print(tsk)
    ttsk = 33.7     # setpoint tsk
    ttcr = 36.5    # setpoint tcr
    ata = 1 # atmospheres?
  
    # Clothing related
    clo = 0.5
    chclo = 1/(0.155*clo)    # pants + t-shirt
    facl = 1.0+0.15*clo     # surface enlargement?

    # Vascular
    skbf = 6.3
    skbfn = 6.3    #forms the equation to calculate new blood flow rate
    skbfl = 90 # L/m2/hr
    multi = 1
    
    # Sweat
    regswl = 500
    csw = 170 # g/m2/hr
    cdil = 75 #litres/(m2/h/K)  #200 in 1986 75 in 1971
    cstr = 0.1                #0.1 in 1986  0.5 in 1971
    sweat = 0

    
    # Body metrics
    ht = phys$ht
    wt = phys$wt
    adu = 0.007184*((ht*100)^0.725)*wt^0.425
    alpha = 0.044 + 0.35/ (skbf-0.1386)
    
    # Metabolism and activity
    mets = 4.5
    rm = ((0.064*wt + 2.896)*11.57)/adu # BMR is for person. equations here use sqm^-1 
#     rm = (63*wt+2896)*0.01157/adu     # kJ to watt schofield
#     rm = (15.057*wt+692.2)*0.04843/adu  # kcal to watt schofield
    #     print (rm)
#     mets = phys$mets
#     rm = 58.2  
    # body fat percentage
#     bmi = wt/ht^2
#     bfp = 1.2*bmi + 0.23*25 - 10.8 - 5.4 #males
#     lbm = wt*(1-(bfp/100)) #lean body mass
# #     katch mcardle formula
#     rm = (370+21.6*lbm)*0.0484/adu
    
    #mifflin st-jour
#     rm = (10*wt+6.25*ht*100-5*25+5)*0.0484/adu
    #harris-benedict
#     rm = (13.397*wt + 4.799*ht*100 - 5.677*25 + 88.362) * 0.0484/adu
#     print(rm)
    
    me = 0.1
    
#    http://deepblue.lib.umich.edu/handle/2027.42/24523
    esk = 7.3 #init
    bz = 0.1

    # Time
    time = 0.0      # init
    exp_time = 30.0 # mins
    hours = exp_time/60
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
    
    
    #takada experiment
#     clo = 0.1 #takada
#     chclo = 1/(0.155*clo)    # pants + t-shirt
#     facl = 1.0+0.15*clo     # surface enlargement?
#     rm = 58.2
#     me = 0
#     mets = 1
#     time = 0.0      # init
#     exp_time = 120.0 # mins
#     hours = exp_time/60
#     interval = 10   # seconds
#     dt = 6
#     steps = exp_time*60/interval
#     data = list()
#     data$TP =  c(rep(29.4,30*dt+1),rep(20,20*dt),rep(29.4,30*dt),rep(40.9,20*dt),rep(29.4,20*dt))  
#     data$MRT =  c(rep(29.6,30*dt+1),rep(20.2,20*dt),rep(29.5,30*dt),rep(40.2,20*dt),rep(29.6,20*dt))
#     data$WS =  c(rep(0.1,30*dt+1),rep(0.24,20*dt),rep(0.10,30*dt),rep(0.12,20*dt),rep(0.12,20*dt))
#     data$RH = c(rep(47.5,30*dt+1),rep(55.6,20*dt),rep(47.6,30*dt),rep(53.5,20*dt),rep(47.8,20*dt))
    

    index = 0
    windex = 1
    while (time < 1.0){
      
#         if (index > 0 && index%%20 == 0){
#             windex = windex+1     #stretch the kestrel data out
#         }
       
      
        if (time >= 0.5) {
          act = rm * mets #no mets     
          me = 0
        }
        else {
          act = rm * mets
          
        }
        index=index+1
        windex=index
        #dynamic weather
        ta = data$TP[windex]
        tr = data$MRT[windex]
        v = data$WS[windex]
        rh = data$RH[windex]/100
        pa = rh * svp(ta)
        #####
        # Radiative and convective
        ####
        
        #still air
        ti = 0.5
        
        #seppanen model
#         if (v < 0.15 && time > 0.5){
#           v=4/3.6
#         }
#         
#         chc = 14.8*v^0.69
      
        
        chc = 4.0*v + 0.35*v * ti - 0.00080*(v*ti)**2 + 3.4 #Ooka
#         chc = 8.4*v^0.5
        chr = 4*0.72*5.67*10^-8*((tcl+tr)/2+273.15)^3
#         chc = 3.1  #takada
#         chr = 4.65 #takada
#         tcl = (chclo*tsk+facl*(chc*ta+chr*tr))/(chclo+facl*(chc+chr))
        
#         while (abs(tcl-tclold) > 0.01){
#             tclold=tcl
#             chr = 4*0.72*5.67*10^-8*((tcl+tr)/2+273.15)^3
#             tcl = (chclo*tsk+facl*(chc*ta+chr*tr))/(chclo+facl*(chc+chr))           
#         }        
        tcl = tsk

        # heat flow from clothing to environment
        rad = facl*chr*(tcl-tr)
#         rad = facl * 0.72 * 0.97 * 5.67*10^-8 * ((tsk+273.15)**4 - (tr+273.15)**4)
        conv = facl*chc*(tcl-ta)
        dry = rad+conv   # R + C negative = gain
#         dry = facl*(chc*(tsk-ta)+chr*(tsk-tr))   # R + C negative = gain

        # dry and latent respiratory heat losses
        eres = 0.017251 * act * (5.8662-pa)      # evaporative loss from respiration
#         eres=0.0023*rm*(44-rh*svp(ta))  #1971 formulation
        
        cres=0.0014 * act * (34.0-ta) # *ata*ff  # convective loss from respiration
        
        wk = act * me
        # heat flow from skin
        hfsk=(tcr-tsk)*(5.28+1.163*skbf)-dry-esk
        
        # heat flow from core
        hfcr=act-(tcr-tsk)*(5.28+1.163*skbf)-cres-eres-wk
        
        # thermal capacities W.HR/C!
        tccr = (1-alpha)*wt*0.97     # 58.2 = 0.97 * 60s!!!!
        tcsk = alpha*wt*0.97
    
        # temperature change per hr (C/hr)
        
        dtsk = (hfsk*adu)/tcsk    
        dtcr = (hfcr*adu)/tccr 
        dtbm=alpha*dtsk+(1.-alpha)*dtcr
        
        dtim = interval/(exp_time*60) # 1 minute / 120 mins
        time = time + dtim
        u = abs(dtsk)
        if(u*dtim > 0.1) {
          dtim = 0.1/u
        }
        u = abs(dtcr)
        if(u*dtim > 0.1) {
          dtim = 0.1/u
        }
        
        tsk = tsk+dtsk*dtim*hours*multi # divided by two because 1.0 time is only 30 mins 
        tcr = tcr+dtcr*dtim*hours*multi

#         print(dtim)
#         print(index)
       
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
        
        skbf = (skbfn+dilat)/(1+stric)
        # limits
        if (skbf < 0.5) {skbf = 0.5}
        if (skbf > skbfl) {skbf = skbfl}

        # skin-core proportion changes with blood flow
        alpha = 0.0417737+ 0.7451832/(skbf+0.58517)
#         alpha = 0.044 + 0.35/ (skbf-0.1386)
#         alpha = 0.1
       

#         cat('warms',warms,'\n')
        #sweating regulation
        regsw = csw*warmb*exp(warms/10.7)
#         regsw = regsw=250.0*warmc+100.0*warmc*warms 
        if (regsw > regswl) {
          regsw = regswl
        }
        
        ersw = 0.68 * regsw
#         ersw = 0.675 * regsw * 2^((tsk-ttsk)/3)   #1971
#         swf = 1.0 #male sweat factor
#         sw = max(8.47 * 10**-5 * ((alpha* tsk + (1-alpha) * tcr) - ttcr),0)*swf
#         sw_h = sw * 3600
#         ersw = sw_h * 675
        
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
        fpcl=1./(1.0+0.143*(chc)*clo)
        
        emax = 2.2*chc*svp_m(tsk)-rh*svp_m(ta)*fpcl
        
        prsw = ersw / emax

        pdif = (1-prsw)*0.06
#         edif = pdif*emax
        vp_m = svp_m(ta)*rh
        edif=-1.694 * 10**-7  * 2430000 * (vp_m - svp_m(tsk)) 
        esk = ersw + edif
#         cat(ersw/0.7,regsw,"\n")
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
        
        sweat = sweat + ((cres+esk+edrip)*adu/0.68)*(interval/3600) #g/timestep
        
        #vapour pressure at skin
        vpsk = pwet * svp(tsk)+(1-pwet)*pa
        # rh at skin
        rhsk = vpsk/svp(tsk)
#         cat("tsk:",tsk,"\n")
        vtsk <- c(vtsk,tsk)
        vdtsk <- c(vdtsk,dtsk)
        vtcr <- c(vtcr,tcr)
#         vtcl <- c(vtcl,tcl)
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
        plot(vtsk,type="l",main="tsk", ylim=c(29,38))#ylim=c(30,36), main= "tsk and tcl (red)")
        lines(vtcr,col="red")
        lines(vtcl,col="blue")
        plot(vtcr,type="l")#,ylim=c(36.5,37.5))
        plot(vr,type="l", main="R")
        plot(vc,type="l", main="C")
        plot(ve,type="l", main = "E",ylim=c(0,120))
        lines(vedif,col="blue")
        lines(versw,col="red")
        
        cat(sweat,'\n')#,"; bmr:",rm,"bfp:",bfp,"lbm:",lbm,"\n")
        return(sweat)
    }

#subjects data
hs1 = list(); hs1$ht = 1.77; hs1$wt = 70.65; hs1$mets1 = 5.56; hs1$mets2 = 5.48;hs1$ttcr1 = 36.73; hs1$ttcr2 = 36.52
hs2 = list(); hs2$ht = 1.78; hs2$wt = 67.765; hs2$mets1 = 3.64; hs2$mets2 = 3.56; hs2$ttcr1 = 36.93; hs2$ttcr2 = 36.87
hs3 = list(); hs3$ht = 1.74; hs3$wt = 70.395; hs3$mets1 = 4.8; hs3$mets2 = 4.84; hs3$ttcr1 = 37.6; hs3$ttcr2 = 36.68
hs4 = list(); hs4$ht = 1.68; hs4$wt = 57.54; hs4$mets1 = 5.8; hs4$mets2 = 5.68; hs4$ttcr1 = 36.37; hs4$ttcr2 = 36.5

ks1 = list(); ks1$ht = 1.8; ks1$wt = 86.515; ks1$mets1 = 3.76; ks1$mets2 = 3.88; ks1$ttcr1 = 36.93; ks1$ttcr2 = 36.93
ks2 = list(); ks2$ht = 1.74; ks2$wt = 72.1; ks2$mets1 = 4; ks2$mets2 = 3.8;ks2$ttcr1 = 36.14; ks2$ttcr2 = 36.54
ks3 = list(); ks3$ht = 1.69; ks3$wt = 53.755; ks3$mets1 = 4.84; ks3$mets2 = 4.64;ks3$ttcr1 = 36.49; ks3$ttcr2 = 36.49
ks4 = list(); ks4$ht = 1.54; ks4$wt = 44.165; ks4$mets1 = 4.88; ks4$mets2 = 4.72;ks4$ttcr1 = 36.64; ks4$ttcr2 = 36.86

ns1 = list(); ns1$ht = 1.69; ns1$wt = 52.705; ns1$mets1 = 5.12; ns1$mets2 = 5.4; ns1$ttcr1 = 36.5;ns1$ttcr2 = 36.5
ns2 = list(); ns2$ht = 1.7; ns2$wt = 71.175; ns2$mets1 = 4.52; ns2$mets2 = 4.08; ns2$ttcr1 = 36.73;ns2$ttcr2 = 36.73
ns3 = list(); ns3$ht = 1.8; ns3$wt = 71.51; ns3$mets1 = 4.36; ns3$mets2 = 4.28; ns3$ttcr1 = 36.5;ns3$ttcr2 = 36.67
ns4 = list(); ns4$ht = 1.69; ns4$wt = 59.98; ns4$mets1 = 3.68; ns4$mets2 = 3.68;ns4$ttcr1 = 37.02;ns4$ttcr2 = 36.62


wrap <- function(exp,candidate,pred){
#   cat(exp,"1\n")
  
  for (n in 1:4){
    phys = list()
    subject = get(paste(candidate,n,sep="")) #ks,hs or ns
    phys$ht = subject$ht
    phys$wt = subject$wt
    phys$mets = subject$mets1
    phys$ttcr = subject$ttcr2
    data<-s2nk("K01",exp,1,30)
    sub(data,phys,exp,n,1) -> x
    pred = c(pred,x)
  }
#   cat("----------------\n")
#   cat(exp,"2\n")
  for (n in 1:4){
    phys = list()
    subject = get(paste(candidate,n,sep="")) #ks,hs or ns
    phys$ht = subject$ht
    phys$wt = subject$wt
    phys$mets = subject$mets2
    phys$ttcr = subject$ttcr2
    data<-s2nk("K01",exp,2,30)
    sub(data,phys,exp,n,2) -> x
    pred = c(pred,x)
  }
  return(pred)
#   cat("----------------\n")
}
pred = c()
wrap("H","hs",pred) -> pred
wrap("K","ks",pred) -> pred
wrap("N","ns",pred) -> pred
# print(pred)
print(summary(lm(pred~sweat$Weight.loss)))
print(rmse(pred,sweat$Weight.loss))

