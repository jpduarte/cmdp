Title "Untitled"

Controls {
}

Definitions {
	Constant "Const.SubstrateFIN" {
		Species = "BoronActiveConcentration"
		Value = 6e+18
	}
	Constant "Const.SubstratePT" {
		Species = "BoronActiveConcentration"
		Value = 2e+18
	}
	Constant "Const.Sext" {
		Species = "ArsenicActiveConcentration"
		Value = 1e+21
	}
	Constant "Const.Dext" {
		Species = "ArsenicActiveConcentration"
		Value = 1e+21
	}
	Refinement "RefDef.Global" {
		MaxElementSize = ( 0.0038 0.0038 0.01 )
		MinElementSize = ( 0.00076 0.00076 0.005 )
	}
}

Placements {
	Constant "PlaceCD.SubstrateFIN" {
		Reference = "Const.SubstrateFIN"
		EvaluateWindow {
			Element = region ["SiFin"]
		}
	}
	Constant "PlaceCD.SubstratePT" {
		Reference = "Const.SubstratePT"
		EvaluateWindow {
			Element = region ["SiSTI"]
		}
	}
	Constant "PlaceCD.Sext" {
		Reference = "Const.Sext"
		EvaluateWindow {
			Element = region ["SiS"]
		}
	}
	Constant "PlaceCD.Dext" {
		Reference = "Const.Dext"
		EvaluateWindow {
			Element = region ["SiD"]
		}
	}
	Refinement "Place.Global" {
		Reference = "RefDef.Global"
		RefineWindow = Cuboid [(-0.015 -0.03 -0.06) (0.015 0.1 0.06)]
	}
}

