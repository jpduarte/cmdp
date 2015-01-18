 
File {
   * input files:
   Grid=   "./v7/finfet_7Lg_5_msh.tdr"
  Parameter = "./v7/models.par"
   * output files:
   Plot=   "@tdrdat@"
   Current="./v7/finfet_7LgWt0.0076Lg1.0Hfin0.042tox0.0008WbR0.0038WbL0.0038vd0.43"
   Output= "@log@"
}

Electrode {
   { Name="source"    Voltage= 0.0 DistResist = 1e-12}
   { Name="drain"     Voltage= 0.0 DistResist = 1e-12 }
   { Name="gate"      Voltage= 0.0 Workfunction= 4.235}
   { Name="substrate" Voltage= 0.0 DistResist = 1e-12}
}



Physics{
	Temperature=300
}

Physics(Material="Silicon"){
   eMobility(
      PhuMob
	  HighFieldSaturation( ParameterSetName = "myset" )
      Enormal
   )
}


Math{
-CheckUndefinedModels
   Method=ILS
   Extrapolate
   Iterations=20
   Notdamped =20
   Number_of_Threads = 4    * For more than 1 CPUs
}


Solve {
  #-initial solution:
  Coupled(Iterations= 100 LineSearchDamping= 1e-4){ Poisson } 
  Coupled { Poisson Electron }
  
  #-- Ramp drain to VdLin
  Quasistationary( 
  InitialStep= 1e-3 Increment= 1.1
    MinStep= 1e-6 MaxStep= 0.01 
    Goal { Name= "drain" Voltage= 0.43 } 
  ){ Coupled {Poisson Electron } }

  #-- Ramp drain to VdLin
  Quasistationary( 
  InitialStep= 1e-3 Increment= 1.1
    MinStep= 1e-6 MaxStep= 0.01 
    Goal { Name= "substrate" Voltage= 0.0 } 
  ){ Coupled {Poisson Electron } }

  #-- Ramp drain to VdLin
  Quasistationary( 
  InitialStep= 1e-3 Increment= 1.1
    MinStep= 1e-6 MaxStep= 0.01 
    Goal { Name= "gate" Voltage= 0.00 } 
  ){ Coupled {Poisson Electron } }

  #-ramp gate:
  Quasistationary ( MaxStep=0.01
                    Goal { Name="gate" Voltage=0.86 } )
                  { Coupled { Poisson Electron } }
}





