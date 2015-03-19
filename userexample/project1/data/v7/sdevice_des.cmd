 
File {
   * input files:
   Grid=   "@tdr@"
  Parameter = "@parameter@"
   * output files:
   Plot=   "@tdrdat@"
   Current="@plot@"
   Output= "@log@"
}

Electrode {
   { Name="source"    Voltage= 0.0 DistResist = 1e-12}
   { Name="drain"     Voltage= 0.0 DistResist = 1e-12 }
   { Name="gate"      Voltage= 0.0 Workfunction= @WF_gate@}
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
    Goal { Name= "drain" Voltage= @Vd_user@ } 
  ){ Coupled {Poisson Electron } }

  #-- Ramp drain to VdLin
  Quasistationary( 
  InitialStep= 1e-3 Increment= 1.1
    MinStep= 1e-6 MaxStep= 0.01 
    Goal { Name= "substrate" Voltage= @Vb_user@ } 
  ){ Coupled {Poisson Electron } }

  #-- Ramp drain to VdLin
  Quasistationary( 
  InitialStep= 1e-3 Increment= 1.1
    MinStep= 1e-6 MaxStep= 0.01 
    Goal { Name= "gate" Voltage= @Vginitial_user@ } 
  ){ Coupled {Poisson Electron } }

  #-ramp gate:
  Quasistationary ( MaxStep=0.01
                    Goal { Name="gate" Voltage=@Vgfinal_user@ } )
                  { Coupled { Poisson Electron } }
}





