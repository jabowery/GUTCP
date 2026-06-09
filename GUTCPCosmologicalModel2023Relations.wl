(* ::Package:: *)

(* GUTCPCosmologicalModel2023Relations.wl

   Relation-first Wolfram Language companion for GUTCPCosmologicalModel.nb.

   The design rule is:

     relation -> constrained fiber -> singleton graph, if forced

   A symbol ending in Rel returns a Boolean relation suitable for Reduce.
   A symbol ending in Expr is algebraic shorthand used inside relations.
   A symbol ending in Reduce asks Mathematica to reduce a relation under stated
   constraints.  ToRules/ReplaceAll should be used only after a branch/domain
   reduction has made the relevant fiber singleton.
*)

ClearAll[
  c, G, sigma, e, mUAtStart, CMBRTempAtStart, tNow, H0Target, CMBTempNow,
  periodNumericExpr, omegaNumericExpr, thetaExpr, radiusMinExpr,
  radiusAmplitudeExpr, radiusExpr, expansionExtentExpr, radiusRateExpr,
  powerExpr, areaExpr, radianceExpr, temperatureExpr, radiusMinFromTempExpr,
  massFromRadiusExpr, massBookExpr, massIntegratedExpr,
  thetaRel, radiusRel, expansionExtentRel, radiusRateRel, hubbleCtRel,
  hubbleRadiusRel, hubbleExpansionExtentRel, powerRel, areaRel, radianceRel,
  temperatureRel, radiusMinFromTempRel, massFromRadiusRel,
  massAtStartFromTempRel, massBookRel, massIntegratedRel,
  calibrationRel, fixedTimeCalibrationRel,
  zCtPathRel, zRadiusScaleRel, zExpansionExtentRel,
  endpointConstantMassRel, endpointBookMassRel, endpointIntegratedMassRel,
  zCombinedExtentEndpointRel, currentExpansionBranchRel,
  radialLightTravelDistanceRel, blueShiftSegmentEntryThresholdRel,
  beyondCurrentExpansionStartRel, symmetricEndpointNeutralThresholdRel,
  symmetricEndpointNetBlueRel,
  currentExpansionDistanceAdmissibleRel, preCurrentSegmentRequiredRel,
  observedRedshiftSignRel, endpointSymmetricSignAdmissibleRel,
  red, neutral, blue,
  cmbrAngularViewExpr, cmbrMultipoleSkyFractionExpr,
  cmbrMultipoleStructureScaleExpr, cmbrMultipoleStructureRel,
  redshiftFiberRel, branchRestrictedFiberRel, emissionFiberReduce,
  singletonFiberRules, thetaExtentRel, thetaExtentCurrentExpansionReduce,
  thetaExtentCurrentExpansionArcSinRel, angularDiameterObsRel,
  nearFlatGeometryClaimRel, tenGyrGeometryCalibrationHazardRel,
  absoluteRestReferenceRel, cepheidQuasiGeometricDistanceRel,
  eModePolarizationRel, bModePolarizationRel, bOverERatioRel,
  gutcpLightBundleConstraintRows,
  gutcpLightBundleRel, angularDiameterGUTCPRel, readinessRows,
  currentExpansion
];

(* ------------------------------------------------------------------------- *)
(* Algebraic expressions from the notebook/numeric period convention.          *)
(* ------------------------------------------------------------------------- *)

periodNumericExpr[m_: mUAtStart] := 2 Pi G m/c^3;
omegaNumericExpr[m_: mUAtStart] := 2 Pi/periodNumericExpr[m];
thetaExpr[t_, m_: mUAtStart] := omegaNumericExpr[m] t;

radiusMinExpr[m_: mUAtStart] := 2 G m/c^2;
radiusAmplitudeExpr[m_: mUAtStart] := 4 Pi G m/c^2;

radiusExpr[t_, m_: mUAtStart] :=
  radiusMinExpr[m] + radiusAmplitudeExpr[m] (1 - Cos[thetaExpr[t, m]]);

expansionExtentExpr[t_, m_: mUAtStart] :=
  radiusExpr[t, m] - radiusMinExpr[m];

radiusRateExpr[t_, m_: mUAtStart] :=
  radiusAmplitudeExpr[m] omegaNumericExpr[m] Sin[thetaExpr[t, m]];

powerExpr[t_, m_: mUAtStart] :=
  c^5 (1 + Cos[thetaExpr[t, m]])/(8 Pi G);

areaExpr[t_, m_: mUAtStart] := 4 Pi radiusExpr[t, m]^2;
radianceExpr[t_, m_: mUAtStart] := powerExpr[t, m]/areaExpr[t, m];
temperatureExpr[t_, m_: mUAtStart] := (radianceExpr[t, m]/(e sigma))^(1/4);

radiusMinFromTempExpr[temp_] :=
  Sqrt[c^5/((4 Pi)^2 G e sigma temp^4)];

massFromRadiusExpr[r_] := c^2 r/(2 G);

massBookExpr[t_, m_: mUAtStart] := m (1 + Cos[thetaExpr[t, m]])/2;

massIntegratedExpr[t_, m_: mUAtStart] :=
  m - c^3 (t + Sin[thetaExpr[t, m]]/omegaNumericExpr[m])/(8 Pi G);

(* ------------------------------------------------------------------------- *)
(* Relation predicates.                                                       *)
(* ------------------------------------------------------------------------- *)

thetaRel[t_, theta_, m_: mUAtStart] := theta == thetaExpr[t, m];
radiusRel[t_, r_, m_: mUAtStart] := r == radiusExpr[t, m];
expansionExtentRel[t_, extent_, m_: mUAtStart] :=
  extent == expansionExtentExpr[t, m];
radiusRateRel[t_, rate_, m_: mUAtStart] := rate == radiusRateExpr[t, m];

(* Mills Eq. 32.156 denominator choice under audit. *)
hubbleCtRel[t_, h_, m_: mUAtStart] := h == radiusRateExpr[t, m]/(c t);

(* Corrected denominator choices. *)
hubbleRadiusRel[t_, h_, m_: mUAtStart] := h == radiusRateExpr[t, m]/radiusExpr[t, m];
hubbleExpansionExtentRel[t_, h_, m_: mUAtStart] :=
  h == radiusRateExpr[t, m]/expansionExtentExpr[t, m];

powerRel[t_, p_, m_: mUAtStart] := p == powerExpr[t, m];
areaRel[t_, area_, m_: mUAtStart] := area == areaExpr[t, m];
radianceRel[t_, rad_, m_: mUAtStart] := rad == radianceExpr[t, m];
temperatureRel[t_, temp_, m_: mUAtStart] := temp == temperatureExpr[t, m];

radiusMinFromTempRel[temp_, rMin_] := rMin == radiusMinFromTempExpr[temp];
massFromRadiusRel[r_, m_] := m == massFromRadiusExpr[r];

massAtStartFromTempRel[temp_, m_] :=
  Exists[{rMin}, radiusMinFromTempRel[temp, rMin] && massFromRadiusRel[rMin, m]];

massBookRel[t_, mass_, m_: mUAtStart] := mass == massBookExpr[t, m];
massIntegratedRel[t_, mass_, m_: mUAtStart] := mass == massIntegratedExpr[t, m];

(* The notebook calibration is a relation on {tObs, tempStart}. *)
calibrationRel[tObs_, tempStart_, h0Target_, cmbTarget_] :=
  Exists[{m, h, tempNow},
    massAtStartFromTempRel[tempStart, m] &&
    hubbleCtRel[tObs, h, m] &&
    temperatureRel[tObs, tempNow, m] &&
    h == h0Target &&
    tempNow == cmbTarget
  ];

(* This relation is intentionally overconstrained if tFixed is fixed and only
   tempStart is left to satisfy both empirical constraints. *)
fixedTimeCalibrationRel[tFixed_, tempStart_, h0Target_, cmbTarget_] :=
  calibrationRel[tFixed, tempStart, h0Target, cmbTarget];

(* ------------------------------------------------------------------------- *)
(* Redshift relations.                                                        *)
(* ------------------------------------------------------------------------- *)

zCtPathRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z ==
    Exp[
      radiusAmplitudeExpr[m] omegaNumericExpr[m]/c
        (SinIntegral[omegaNumericExpr[m] tObs] -
         SinIntegral[omegaNumericExpr[m] tEmit])
    ];

zRadiusScaleRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z == radiusExpr[tObs, m]/radiusExpr[tEmit, m];

zExpansionExtentRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z == expansionExtentExpr[tObs, m]/expansionExtentExpr[tEmit, m];

endpointConstantMassRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z ==
    (1 + 2 G m/(c^2 radiusExpr[tEmit, m]))/
    (1 + 2 G m/(c^2 radiusExpr[tObs, m]));

endpointBookMassRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z ==
    (1 + 2 G massBookExpr[tEmit, m]/(c^2 radiusExpr[tEmit, m]))/
    (1 + 2 G massBookExpr[tObs, m]/(c^2 radiusExpr[tObs, m]));

endpointIntegratedMassRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z ==
    (1 + 2 G massIntegratedExpr[tEmit, m]/(c^2 radiusExpr[tEmit, m]))/
    (1 + 2 G massIntegratedExpr[tObs, m]/(c^2 radiusExpr[tObs, m]));

zCombinedExtentEndpointRel[tEmit_, tObs_, z_, m_: mUAtStart] :=
  1 + z ==
    (expansionExtentExpr[tObs, m]/expansionExtentExpr[tEmit, m]) *
    (1 + 2 G m/(c^2 radiusExpr[tEmit, m]))/
    (1 + 2 G m/(c^2 radiusExpr[tObs, m]));

(* Branch information is explicitly separate from the redshift relation. *)
currentExpansionBranchRel[tEmit_, tObs_, m_: mUAtStart] :=
  0 <= thetaExpr[tEmit, m] <= thetaExpr[tObs, m] < Pi;

(* Radial light-travel distance is itself a relation.  If tObs is measured
   from the start of the current expansion, dLight == c tObs is the threshold
   at which the emission event leaves the current-expansion branch. *)
radialLightTravelDistanceRel[tEmit_, tObs_, dLight_] :=
  dLight == c (tObs - tEmit);

blueShiftSegmentEntryThresholdRel[dThreshold_, tObs_] :=
  dThreshold == c tObs;

beyondCurrentExpansionStartRel[dLight_, tObs_] :=
  dLight > c tObs;

(* For the even oscillatory endpoint relations R(-t)==R(t) and E(-t)==E(t),
   the endpoint-only net-red/net-blue boundary occurs at tEmit == -tObs, i.e.
   dLight == 2 c tObs.  This is distinct from the threshold at which a
   pre-current-expansion blueshift segment first enters the path. *)
symmetricEndpointNeutralThresholdRel[dThreshold_, tObs_] :=
  dThreshold == 2 c tObs;

symmetricEndpointNetBlueRel[dLight_, tObs_] :=
  dLight > 2 c tObs;

currentExpansionDistanceAdmissibleRel[dLight_, tObs_] :=
  0 <= dLight <= c tObs;

preCurrentSegmentRequiredRel[dLight_, tObs_] :=
  dLight > c tObs;

observedRedshiftSignRel[z_, sign_] :=
  (sign == red && z > 0) ||
  (sign == neutral && z == 0) ||
  (sign == blue && -1 < z < 0);

endpointSymmetricSignAdmissibleRel[dLight_, tObs_, sign_] :=
  (sign == red && dLight < 2 c tObs) ||
  (sign == neutral && dLight == 2 c tObs) ||
  (sign == blue && dLight > 2 c tObs);

(* Mills page-1579 CMBR multipole structure-scale relation, not BAO. *)
cmbrAngularViewExpr[rSphereLy_] := 2 Pi rSphereLy;
cmbrMultipoleSkyFractionExpr[ell_] := 2/ell;
cmbrMultipoleStructureScaleExpr[ell_, rSphereLy_] :=
  cmbrAngularViewExpr[rSphereLy] cmbrMultipoleSkyFractionExpr[ell];

cmbrMultipoleStructureRel[ell_, rSphereLy_, scaleLy_] :=
  ell > 0 && scaleLy == cmbrMultipoleStructureScaleExpr[ell, rSphereLy];

redshiftFiberRel[redshiftRel_, tEmit_, tObs_, z_, m_: mUAtStart] :=
  tEmit < tObs && redshiftRel[tEmit, tObs, z, m];

branchRestrictedFiberRel[
  redshiftRel_, branchRel_, tEmit_, tObs_, z_, m_: mUAtStart
] :=
  redshiftFiberRel[redshiftRel, tEmit, tObs, z, m] &&
  branchRel[tEmit, tObs, m];

emissionFiberReduce[
  redshiftRel_, branchRel_, zValue_, tObsValue_, tEmit_Symbol,
  mValue_: mUAtStart, assumptions_: True
] :=
  Reduce[
    branchRestrictedFiberRel[redshiftRel, branchRel, tEmit, tObsValue, zValue, mValue] &&
      assumptions,
    tEmit,
    Reals
  ];

(* Use only after Reduce has established a singleton fiber. *)
singletonFiberRules[reducedRel_] := ToRules[reducedRel];

(* ------------------------------------------------------------------------- *)
(* Arcsine as relation, then branch-restricted singleton representation.       *)
(* ------------------------------------------------------------------------- *)

thetaExtentRel[thetaEmit_, z_, tObs_, m_: mUAtStart] :=
  Sin[thetaEmit/2]^2 ==
    expansionExtentExpr[tObs, m]/((1 + z) 2 radiusAmplitudeExpr[m]);

thetaExtentCurrentExpansionReduce[
  zValue_, tObsValue_, thetaEmit_Symbol, mValue_: mUAtStart,
  assumptions_: True
] :=
  Reduce[
    thetaExtentRel[thetaEmit, zValue, tObsValue, mValue] &&
      0 <= thetaEmit <= thetaExpr[tObsValue, mValue] < Pi &&
      assumptions,
    thetaEmit,
    Reals
  ];

(* This is not the primitive relation.  It is the singleton-branch
   representation after the current-expansion branch constraints are admitted. *)
thetaExtentCurrentExpansionArcSinRel[thetaEmit_, z_, tObs_, m_: mUAtStart] :=
  0 <= thetaEmit <= thetaExpr[tObs, m] < Pi &&
  thetaEmit ==
    2 ArcSin[
      Sqrt[
        expansionExtentExpr[tObs, m]/((1 + z) 2 radiusAmplitudeExpr[m])
      ]
    ];

(* ------------------------------------------------------------------------- *)
(* Angular-diameter relation and the deliberately missing GUTCP optical law.   *)
(* ------------------------------------------------------------------------- *)

angularDiameterObsRel[z_, dA_, sourceArea_, observedSolidAngle_] :=
  dA^2 == sourceArea/observedSolidAngle;

(* Constraints on the still-open GUTCP light-bundle relation from reconstructed
   pages 1569, 1577-1579, and the page-1602 footnote.  These are not a d_A
   closure by themselves. *)
nearFlatGeometryClaimRel[cmbrAngularSizeDeg_, tExpansionYears_] :=
  cmbrAngularSizeDeg == 1 && tExpansionYears == 10^10;

tenGyrGeometryCalibrationHazardRel[tClaimYears_, tCalibratedYears_] :=
  tClaimYears == 10^10 && tCalibratedYears != tClaimYears;

absoluteRestReferenceRel[translationalVelocity_, correctedVelocity_] :=
  translationalVelocity == 0 || correctedVelocity == 0;

cepheidQuasiGeometricDistanceRel[
  distance_, angularDiameterChange_, radialVelocityCurve_, pulsationCycle_
] :=
  distance > 0 &&
  angularDiameterChange > 0 &&
  radialVelocityCurve =!= Missing["Unspecified"] &&
  pulsationCycle > 0;

eModePolarizationRel[ell_, deltaTE_, cEff_] :=
  ell > 0 &&
  deltaTE == cEff 77 Sinc[(Pi/140) (ell + 197)];

bModePolarizationRel[ell_, deltaTB_, cEff_, bOverERatio_] :=
  ell > 0 &&
  deltaTB == bOverERatio cEff 77 Sinc[(Pi/140) (ell + 197 + 70)];

bOverERatioRel[bOverERatio_, deltaAleph_, lightSphereRadius_] :=
  bOverERatio == deltaAleph/lightSphereRadius;

gutcpLightBundleConstraintRows = {
  {"p.1569", "nearFlatGeometryClaimRel, tenGyrGeometryCalibrationHazardRel", "near-flat CMBR angular-structure claim tied to fixed 10 Gyr approximation"},
  {"p.1602 footnote", "absoluteRestReferenceRel, cepheidQuasiGeometricDistanceRel", "absolute-rest r-sphere reference and quasi-geometrical Cepheid distance calibration"},
  {"p.1577 Eq.32.203", "eModePolarizationRel", "E-mode CMBR multipole phase/amplitude"},
  {"p.1577 Eqs.32.204-32.205", "bModePolarizationRel, bOverERatioRel", "E-to-B conversion and DeltaAleph/(c t) amplitude ratio"},
  {"p.1578", "lensing/ring-arc prose relation", "photon-bundle lensing to observed angular structures"},
  {"p.1579", "cmbrMultipoleStructureRel", "CMBR multipole ruler relation, not BAO"}
};

(* gutcpLightBundleRel is intentionally left undefined.  It is the relation
   that must be supplied before the Koksbang-Heinesen C/O/M diagnostics can be
   reduced fairly for GUTCP.  The rows above constrain this frontier without
   closing it. *)

angularDiameterGUTCPRel[
  z_, dA_, tEmit_, tObs_, sourceArea_, observedSolidAngle_,
  branchTag_: currentExpansion
] :=
  angularDiameterObsRel[z, dA, sourceArea, observedSolidAngle] &&
  gutcpLightBundleRel[branchTag, tEmit, tObs, sourceArea, observedSolidAngle];

readinessRows = {
  {"z relation", "zExpansionExtentRel, zRadiusScaleRel, zCtPathRel", "stated"},
  {"branch fiber", "branchRestrictedFiberRel", "stated"},
  {"blueshift segment threshold", "blueShiftSegmentEntryThresholdRel", "stated"},
  {"net endpoint-blue threshold", "symmetricEndpointNeutralThresholdRel", "conditional"},
  {"red/blue branch pruning", "endpointSymmetricSignAdmissibleRel", "conditional"},
  {"CMBR multipole structure scale", "cmbrMultipoleStructureRel", "stated; not BAO"},
  {"singleton extraction", "singletonFiberRules", "only after Reduce"},
  {"H_parallel relation", "hubbleExpansionExtentRel or selected denominator", "conditional"},
  {"d_A observation", "angularDiameterObsRel", "stated"},
  {"light-bundle frontier constraints", "gutcpLightBundleConstraintRows", "stated; closure open"},
  {"BAO ruler relation", "D_M/r_d, H_parallel r_d, or D_V/r_d rather than bare d_A", "conditional"},
  {"GUTCP angular law", "gutcpLightBundleRel", "intentionally unspecified"},
  {"C/O/M reduction", "requires d_A relation joined to optical law", "not yet fair"}
};
