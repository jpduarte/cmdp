(sde:clear)
(sdegeo:set-default-boolean "BAB")

(define Lg @Lg_user@)
(define Lg_2_l (* Lg -0.5))
(define Lg_2_r (* Lg 0.5))

(define Lg_ext 0.010)

(define Lg_S_l (- Lg_2_l Lg_ext))
(define Lg_S_r Lg_2_l)

(define Lg_D_l Lg_2_r)
(define Lg_D_r (+ Lg_2_r Lg_ext))

(define Wsub 0.03)
(define Wsub_l (* Wsub -0.5))
(define Wsub_r (* Wsub 0.5))

(define Hgate_t 0.1)
(define Hgate_b 0.0)

(define Hfin_t @Hfin_user@)

(define bottom_structure 0.03)

(define Wsti Wsub)
(define Wsti_l (* Wsti -0.5))
(define Wsti_r (* Wsti 0.5))

(define Hsti_t 0)
(define Hsti_b (* bottom_structure -1))

(define Hsipt_t 0)
(define Hsipt_b (* bottom_structure -1))

(define Ro  @Ro_user@) 
;--------------device size definitions--------------------
(define Hsi Hfin_t)
(define WsiT @Wfintop_user@)
(define WsiBL @WsiBL_user@)
(define WsiBR @WsiBR_user@)

;TODO: check toxR
(define toxT @tox_user@)
(define toxL @tox_user@)
(define toxR @tox_user@)

(define Wsipt_l (* (+ (* WsiT 0.5) WsiBL) -1.0))
(define Wsipt_r (+ (* WsiT 0.5) WsiBR))

;--------------coordinate internal definitions--------------------
(define sin_alpha (/ Hsi (sqrt (+ (expt WsiBL 2) (expt Hsi 2)))))
(define tan_alpha (/ Hsi WsiBL))
(define cos_alpha (/ WsiBL (sqrt (+ (expt WsiBL 2) (expt Hsi 2)))))

(define sin_alpha2 (/ Hsi (sqrt (+ (expt WsiBR 2) (expt Hsi 2)))))
(define tan_alpha2 (/ Hsi WsiBR))
(define cos_alpha2 (/ WsiBR (sqrt (+ (expt WsiBR 2) (expt Hsi 2)))))

(define Pfin1_X (+ (/ WsiT 2 ) WsiBL))
(define Pfin1_Y 0)

(define Pfin2_X (/ WsiT 2 ))
(define Pfin2_Y Hsi)

(define Pfin3_X (+ (/ WsiT 2 ) WsiBR)) 
(define Pfin3_Y 0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define Ptox1_X (+ WsiBL (+ (/ WsiT 2) (/ toxL sin_alpha))))
(define Ptox1_Y 0)

(define Ptox2_X (+ (/ WsiT 2) (/ toxL sin_alpha)))
(define Ptox2_Y Hsi)

(define Ptox3_X (/ WsiT 2))
(define Ptox3_Y (+ Hsi toxL))

(define Ptox4_X (/ WsiT 2))
(define Ptox4_Y (+ Hsi toxL))

(define Ptox5_X (+ (/ WsiT 2) (/ toxL sin_alpha2)))
(define Ptox5_Y Hsi)


(define Ptox6_X (+ WsiBR (+ (/ WsiT 2) (/ toxL sin_alpha2))))
(define Ptox6_Y 0)

;----------------------------------------------------------------------
; Silicon PT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:create-cuboid (position  (* Pfin1_X -1) Hsipt_b Lg_2_l) (position Pfin3_X Hsipt_t Lg_2_r)  
 "Silicon" "SiSTI")
 
; Oxide ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:create-cuboid (position  Wsti_l Hsti_b Lg_2_l) (position Wsti_r Hsti_t Lg_2_r)  
 "Oxide" "OxideSTI")
 
; Si Fin ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:create-polygon
(list (position (* Pfin1_X -1) Pfin1_Y Lg_2_l)
      (position (* Pfin2_X -1) Pfin2_Y Lg_2_l)
      (position Pfin2_X Pfin2_Y Lg_2_l)
      (position Pfin3_X Pfin3_Y Lg_2_l) 
      (position (* Pfin1_X -1) Pfin1_Y Lg_2_l) ) 
      "Silicon" "SiFin")
(sdegeo:extrude (list (car (find-face-id (position 0.00 (/ Hsi 2) Lg_2_l ))))  Lg ) 

; (sdegeo:fillet-edges 
; (list (car (find-edge-id (position (* Pfin2_X -1) Pfin2_Y 0))) 
;      (car (find-edge-id( position Pfin2_X Pfin2_Y 0))) )  Ro)

; Oxide ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:create-polygon
(list (position (* Ptox1_X -1) Ptox1_Y (- Lg_2_l 0.001) )
       (position (* Ptox2_X -1) Ptox2_Y (- Lg_2_l 0.001) )
       (position (* Ptox3_X -1) Ptox3_Y (- Lg_2_l 0.001) )
       (position Ptox4_X Ptox4_Y (- Lg_2_l 0.001) )
       (position Ptox5_X Ptox5_Y (- Lg_2_l 0.001) ) 
       (position Ptox6_X Ptox6_Y (- Lg_2_l 0.001) ) 
       (position (* Ptox1_X -1) Ptox1_Y (- Lg_2_l 0.001) )) 
       "Oxide" "OxideGate")
(sdegeo:extrude (list (car (find-face-id (position 0.00 (+ Hsi (/ toxT 2)) (- Lg_2_l 0.001)  ))))  (+ Lg 0.002) ) 

; (sdegeo:fillet-edges 
; (list (car (find-edge-id (position (* Ptox2_X -1) Ptox2_Y 0))) 
;      (car (find-edge-id( position Ptox2_X Ptox2_Y 0))) )  Ro)

; Gate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:create-cuboid (position  Wsub_l Hgate_b Lg_2_l) (position Wsub_r Hgate_t Lg_2_r) "Metal" "MetalGate")
 
;; SOURCCE and DRAIN;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(sdegeo:set-default-boolean "ABA")
(sdegeo:create-polygon
(list (position (* Pfin1_X -1) Pfin1_Y Lg_S_l)
      (position (* Pfin2_X -1) Pfin2_Y Lg_S_l)
      (position Pfin2_X Pfin2_Y Lg_S_l)
      (position Pfin3_X Pfin3_Y Lg_S_l) 
      (position (* Pfin1_X -1) Pfin1_Y Lg_S_l) ) 
      "Silicon" "SiS")
(sdegeo:extrude (list (car (find-face-id (position 0.00 (/ Hsi 2) Lg_S_l ))))  Lg_ext ) 
  
;(sdegeo:fillet-edges 
;(list (car (find-edge-id (position (* Pfin2_X -1) Pfin2_Y (+ Lg_S_l 0.001)))) 
;      (car (find-edge-id( position Pfin2_X Pfin2_Y (+ Lg_S_l 0.001) ))) )  Ro)
  
(sdegeo:create-polygon
(list (position (* Pfin1_X -1) Pfin1_Y Lg_D_r)
      (position (* Pfin2_X -1) Pfin2_Y Lg_D_r)
      (position Pfin2_X Pfin2_Y Lg_D_r)
      (position Pfin3_X Pfin3_Y Lg_D_r) 
      (position (* Pfin1_X -1) Pfin1_Y Lg_D_r) ) 
     "Silicon" "SiD")
(sdegeo:extrude (list (car (find-face-id (position  0.00 (/ Hsi 2) Lg_D_r ))))  (* Lg_ext -1) ) 
       
;(sdegeo:fillet-edges 
;(list (car (find-edge-id (position (* Pfin2_X -1) Pfin2_Y (- Lg_D_r 0.001)))) 
;       (car (find-edge-id( position Pfin2_X Pfin2_Y (- Lg_D_r 0.001) ))) )  Ro)
;;;;;;;;;;;;;;;;;;;;;;;;
(sdedr:define-constant-profile "Const.SubstrateFIN" 
 "BoronActiveConcentration" @Doping_FIN@ )
(sdedr:define-constant-profile-region  "PlaceCD.SubstrateFIN" 
 "Const.SubstrateFIN" "SiFin" )
 
(sdedr:define-constant-profile "Const.SubstratePT" "BoronActiveConcentration" 2e18)
(sdedr:define-constant-profile-region  "PlaceCD.SubstratePT" "Const.SubstratePT" "SiSTI" )
 
(sdedr:define-constant-profile "Const.Sext" 
 "ArsenicActiveConcentration" @Doping_PUNCH@ )
(sdedr:define-constant-profile-region  "PlaceCD.Sext" 
"Const.Sext" "SiS" )
 
(sdedr:define-constant-profile "Const.Dext" 
 "ArsenicActiveConcentration" @DOPING_SD@ )
(sdedr:define-constant-profile-region  "PlaceCD.Dext" 
 "Const.Dext" "SiD" )
 ;;;;;;;;;;;;;
 
 ;------------------Define contacts------------------------------------------------------------------------------
(sdegeo:define-contact-set "source" 4.0 (color:rgb 1.0 0.0 0.0) "##")
(sdegeo:define-contact-set "drain" 4.0 (color:rgb 0.0 1.0 0.0) "==")
(sdegeo:define-contact-set "gate" 4.0 (color:rgb 0.0 0.0 1.0) "||")
(sdegeo:define-contact-set "substrate" 4.0 (color:rgb 1.0 0.0 1.0) "<><>")
;----------GATE-------------------------------------------------------------------------------------------------

;OPTION 2 for GATE
(sdegeo:set-current-contact-set "gate")
(define GateID (find-body-id (position 0 Hgate_t 0)))
(sdegeo:set-contact-boundary-faces GateID)
(sdegeo:delete-region GateID)

;----------SOURCE-----------------------------------------------------------------------------------------------
(sdegeo:set-current-contact-set "source")
(sdegeo:set-contact-faces (find-face-id (position 0.0 (/ Hfin_t 2) Lg_S_l)))

;----------DRAIN------------------------------------------------------------------------------------------------
(sdegeo:set-current-contact-set "drain")
(sdegeo:set-contact-faces (find-face-id (position 0.0 (/ Hfin_t 2) Lg_D_r)))
;----------substrate------------------------------------------------------------------------------------------------
(sdegeo:set-current-contact-set "substrate")
(sdegeo:set-contact-faces (find-face-id (position  0.00 Hsipt_b 0 )))
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;MESH;;;;;;;;;;;;;;;;;
; Meshing strategies
(sdedr:define-refinement-window "RefWin.Global" "Cuboid"
(position  Wsub_l Hsti_b Lg_S_l) (position Wsub_r Hgate_t Lg_D_r)  )

(sdedr:define-refinement-size "RefDef.Global"
  (* WsiT 0.5)  (* WsiT 0.5)  (* 0.1 Lg)
  (* WsiT 0.1) (* WsiT 0.1) (* 0.05 Lg))

(sdedr:define-refinement-placement "Place.Global"
 "RefDef.Global" "RefWin.Global" )
 
 ; Multiboxes
;(sdedr:define-refinement-window "RefWin.Channel"
; "Cuboid" (position (* Ptox1_X -1.0) Ptox1_Y Lg_S_l) (position Ptox1_X Ptox2_Y Lg_D_r))
 
;(sdedr:define-multibox-size "RefDefMB.Channel"
; 0.001  0.001  0.01
; 0.0001 0.0001 0.001
 ;1    1   1 )


; Meshing structure 
 ;----------------------------------------------------------------------
; Save CMD file
(sde:build-mesh "snmesh" " " "n@node@_msh")
