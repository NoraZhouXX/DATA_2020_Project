%path = @runpath
cd %path
load raw_data.wf1

pageselect data

smpl @all

genr log_GDP       = @log(rgdp)*100
genr log_cpi         =  @log(pce)*100 
genr shock          = ffr

genr log_connectedness  = @log(c_index2)*100
log_connectedness.hpf c_trend @ c_cycle

series pshock  = @recode(c_cycle >=0,  1, 0)
series nshock  = @recode(c_cycle <0,   1, 0)

series p1  = pshock(-1)
series n1  = nshock(-1)


table est_1
est_1(1,1) = "point_linear"
est_1(1,2) = "point_high"
est_1(1,3) = "point_low"
est_1(1,4) = "linear1"
est_1(1,5) = "linear2"
est_1(1,6) = "linear3"
est_1(1,7) = "H1"
est_1(1,8) = "H2"
est_1(1,9) = "H3"
est_1(1,10) = "L1"
est_1(1,11) = "L2"
est_1(1,12) = "L3"

table est_2
est_2(1,1) = "point_linear"
est_2(1,2) = "point_high"
est_2(1,3) = "point_low"
est_2(1,4) = "linear1"
est_2(1,5) = "linear2"
est_2(1,6) = "linear3"
est_2(1,7) = "H1"
est_2(1,8) = "H2"
est_2(1,9) = "H3"
est_2(1,10) = "L1"
est_2(1,11) = "L2"
est_2(1,12) = "L3"

table est_3
est_3(1,1) = "point_linear"
est_3(1,2) = "point_high"
est_3(1,3) = "point_low"
est_3(1,4) = "linear1"
est_3(1,5) = "linear2"
est_3(1,6) = "linear3"
est_3(1,7) = "H1"
est_3(1,8) = "H2"
est_3(1,9) = "H3"
est_3(1,10) = "L1"
est_3(1,11) = "L2"
est_3(1,12) = "L3"

scalar conf = 1  

smpl 1981q1 2007q4

for !h=0 to 20
	equation lpm1.ls(cov=hac) log_gdp(+!h)       c  ffr    ffr(-1 to -4)  log_GDP(0 to -4)   log_cpi(0 to -4)  @trend @trend^2                             
	equation lpm2.ls(cov=hac) log_cpi(+!h)        c   ffr    ffr(-1 to -4)  log_GDP(0 to -4)   log_cpi(0 to -4)  @trend @trend^2      
	equation lpm3.ls(cov=hac) ffr(+!h)                c   ffr    ffr(-1 to -4)  log_GDP(0 to -4)   log_cpi(0 to -4)   @trend @trend^2                

	equation adlpm1.ls(cov=hac)  log_gdp(+!h)  p1*shock n1*shock  _
	 p1    log_GDP(0)*p1 log_GDP(-1)*p1  log_GDP(-2)*p1  log_GDP(-3)*p1  log_GDP(-4)*p1  ffr(-1)*p1  ffr(-2)*p1  ffr(-3)*p1  ffr(-4)*p1  log_cpi(0)*p1 log_cpi(-1)*p1  log_cpi(-2)*p1  log_cpi(-3)*p1  log_cpi(-4)*p1  _
 	 n1     log_GDP(0)*n1 log_GDP(-1)*n1  log_GDP(-2)*n1  log_GDP(-3)*n1  log_GDP(-4)*n1  ffr(-1)*n1  ffr(-2)*n1  ffr(-3)*n1  ffr(-4)*n1  log_cpi(0)*n1 log_cpi(-1)*n1  log_cpi(-2)*n1  log_cpi(-3)*n1  log_cpi(-4)*n1  @trend @trend^2                                                          

	equation adlpm2.ls(cov=hac)  log_cpi(+!h)  p1*shock n1*shock  _
	 p1  log_GDP(0)*p1  log_GDP(-1)*p1  log_GDP(-2)*p1  log_GDP(-3)*p1  log_GDP(-4)*p1  ffr(-1)*p1  ffr(-2)*p1  ffr(-3)*p1  ffr(-4)*p1   log_cpi(0)*p1 log_cpi(-1)*p1  log_cpi(-2)*p1  log_cpi(-3)*p1  log_cpi(-4)*p1  _
 	 n1  log_GDP(0)*n1 log_GDP(-1)*n1  log_GDP(-2)*n1  log_GDP(-3)*n1  log_GDP(-4)*n1  ffr(-1)*n1  ffr(-2)*n1  ffr(-3)*n1  ffr(-4)*n1   log_cpi(0)*n1  log_cpi(-1)*n1  log_cpi(-2)*n1  log_cpi(-3)*n1  log_cpi(-4)*n1   @trend @trend^2         

	equation adlpm3.ls(cov=hac)  ffr(+!h)  p1*shock n1*shock  _
	 p1   ffr(-1)*p1   ffr(-2)*p1   ffr(-3)*p1   ffr(-4)*p1 log_GDP(0)*p1  log_GDP(-1)*p1  log_GDP(-2)*p1  log_GDP(-3)*p1  log_GDP(-4)*p1 log_cpi(0)*p1  log_cpi(-1)*p1  log_cpi(-2)*p1  log_cpi(-3)*p1  log_cpi(-4)*p1 _
 	 n1      ffr(-1)*n1   ffr(-2)*n1   ffr(-3)*n1   ffr(-4)*n1  log_GDP(0)*n1  log_GDP(-1)*n1  log_GDP(-2)*n1  log_GDP(-3)*n1  log_GDP(-4)*n1  log_cpi(0)*n1 log_cpi(-1)*n1  log_cpi(-2)*n1  log_cpi(-3)*n1  log_cpi(-4)*n1  @trend @trend^2  


est_1(!h+2,1)  =  -1*(lpm1.@coef(2))       
est_1(!h+2,2)  =  -1*adlpm1.@coef(1)   
est_1(!h+2,3)  =  -1*adlpm1.@coef(2) 
est_1(!h+2,4)  =  -1*(lpm1.@coef(2))   - conf*lpm1.@stderrs(2) 
est_1(!h+2,5)  =  -1*(lpm1.@coef(2))   
est_1(!h+2,6)  =  -1*(lpm1.@coef(2))   + conf*lpm1.@stderrs(2)  
est_1(!h+2,7)  = -1*adlpm1.@coef(1) - conf*adlpm1.@stderrs(1)  
est_1(!h+2,8)  =  -1*adlpm1.@coef(1)      
est_1(!h+2,9)  =   -1*adlpm1.@coef(1) + conf*adlpm1.@stderrs(1)  
est_1(!h+2,10) =  -1*adlpm1.@coef(2) - conf*adlpm1.@stderrs(2) 
est_1(!h+2,11) =    -1*adlpm1.@coef(2)    
est_1(!h+2,12) =    -1*adlpm1.@coef(2) + conf* adlpm1.@stderrs(2)  


est_2(!h+2,1)  =  -1*(lpm2.@coef(2))       
est_2(!h+2,2)  =  -1*adlpm2.@coef(1)   
est_2(!h+2,3)  =  -1*adlpm2.@coef(2) 
est_2(!h+2,4)  =  -1*(lpm2.@coef(2))   - conf*lpm2.@stderrs(2) 
est_2(!h+2,5)  =  -1*(lpm2.@coef(2))   
est_2(!h+2,6)  =  -1*(lpm2.@coef(2))   + conf*lpm2.@stderrs(2)  
est_2(!h+2,7)  = -1*adlpm2.@coef(1) - conf*adlpm2.@stderrs(1)  
est_2(!h+2,8)  =  -1*adlpm2.@coef(1)      
est_2(!h+2,9)  =   -1*adlpm2.@coef(1) + conf*adlpm2.@stderrs(1)  
est_2(!h+2,10) =  -1*adlpm2.@coef(2) - conf*adlpm2.@stderrs(2) 
est_2(!h+2,11) =    -1*adlpm2.@coef(2)    
est_2(!h+2,12) =    -1*adlpm2.@coef(2) + conf* adlpm2.@stderrs(2)  


est_3(!h+2,1)  =  -1*(lpm3.@coef(2))       
est_3(!h+2,2)  =  -1*adlpm3.@coef(1)   
est_3(!h+2,3)  =  -1*adlpm3.@coef(2) 
est_3(!h+2,4)  =  -1*(lpm3.@coef(2))   - conf*lpm3.@stderrs(2) 
est_3(!h+2,5)  =  -1*(lpm3.@coef(2))   
est_3(!h+2,6)  =  -1*(lpm3.@coef(2))   + conf*lpm3.@stderrs(2)  
est_3(!h+2,7)  = -1*adlpm3.@coef(1) - conf*adlpm3.@stderrs(1)  
est_3(!h+2,8)  =  -1*adlpm3.@coef(1)      
est_3(!h+2,9)  =   -1*adlpm3.@coef(1) + conf*adlpm3.@stderrs(1)  
est_3(!h+2,10) =  -1*adlpm3.@coef(2) - conf*adlpm3.@stderrs(2) 
est_3(!h+2,11) =    -1*adlpm3.@coef(2)    
est_3(!h+2,12) =    -1*adlpm3.@coef(2) + conf* adlpm3.@stderrs(2)  

next

est_1.save(t=csv) est_1.csv
est_2.save(t=csv) est_2.csv
est_3.save(t=csv) est_3.csv


