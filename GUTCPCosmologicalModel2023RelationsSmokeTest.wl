(* ::Package:: *)

(* Smoke test for GUTCPCosmologicalModel2023Relations.wl.

   This is intended to be evaluated either from the Wolfram front end or from a
   batch kernel when one is available.  It deliberately tests relational shape,
   not numerical calibration.
*)

ClearAll[relationFileDirectory, relationFile, relationSmokeReport];

relationFileDirectory[] := Module[{fromInput, fromNotebook},
  fromInput = Quiet[Check[DirectoryName[$InputFileName], $Failed]];
  If[StringQ[fromInput] && fromInput =!= "", Return[fromInput]];

  fromNotebook = Quiet[Check[NotebookDirectory[], $Failed]];
  If[StringQ[fromNotebook] && fromNotebook =!= "", Return[fromNotebook]];

  Directory[]
];

relationFile =
  FileNameJoin[{relationFileDirectory[], "GUTCPCosmologicalModel2023Relations.wl"}];

Get[relationFile];

relationSmokeReport = <|
  "relationFile" -> relationFile,
  "zExpansionExtentRelIsEquation" ->
    MatchQ[zExpansionExtentRel[tEmit, tObs, z, m], _Equal],
  "zRadiusScaleRelIsEquation" ->
    MatchQ[zRadiusScaleRel[tEmit, tObs, z, m], _Equal],
  "redshiftFiberRelIsConjunction" ->
    MatchQ[redshiftFiberRel[zExpansionExtentRel, tEmit, tObs, z, m], _And],
  "branchRestrictedFiberRelIsConjunction" ->
    MatchQ[
      branchRestrictedFiberRel[
        zExpansionExtentRel,
        currentExpansionBranchRel,
        tEmit,
        tObs,
        z,
        m
      ],
      _And
    ],
  "radialLightTravelDistanceRelIsEquation" ->
    MatchQ[radialLightTravelDistanceRel[tEmit, tObs, dLight], _Equal],
  "blueShiftSegmentEntryThresholdRelIsEquation" ->
    MatchQ[blueShiftSegmentEntryThresholdRel[dThreshold, tObs], _Equal],
  "symmetricEndpointNeutralThresholdRelIsEquation" ->
    MatchQ[symmetricEndpointNeutralThresholdRel[dThreshold, tObs], _Equal],
  "symmetricEndpointNetBlueRelIsInequality" ->
    MatchQ[symmetricEndpointNetBlueRel[dLight, tObs], _Greater],
  "currentExpansionDistanceAdmissibleRelIsInequality" ->
    MatchQ[currentExpansionDistanceAdmissibleRel[dLight, tObs], _Inequality],
  "preCurrentSegmentRequiredRelIsInequality" ->
    MatchQ[preCurrentSegmentRequiredRel[dLight, tObs], _Greater],
  "observedRedshiftSignRelIsDisjunction" ->
    MatchQ[observedRedshiftSignRel[z, sign], _Or],
  "endpointSymmetricBlueSignAdmissibleRelIsInequality" ->
    MatchQ[endpointSymmetricSignAdmissibleRel[dLight, tObs, blue], _Greater],
  "cmbrMultipoleStructureRelIsConjunction" ->
    MatchQ[cmbrMultipoleStructureRel[ell, rSphereLy, scaleLy], _And],
  "thetaExtentRelIsEquation" ->
    MatchQ[thetaExtentRel[thetaEmit, z, tObs, m], _Equal],
  "arcSinRepresentationIsBranchConditional" ->
    MatchQ[
      thetaExtentCurrentExpansionArcSinRel[thetaEmit, z, tObs, m],
      _And
    ],
  "angularDiameterObsRelIsEquation" ->
    MatchQ[angularDiameterObsRel[z, dA, sourceArea, observedSolidAngle], _Equal],
  "nearFlatGeometryClaimRelIsConjunction" ->
    MatchQ[nearFlatGeometryClaimRel[cmbrAngularSizeDeg, tExpansionYears], _And],
  "tenGyrGeometryCalibrationHazardRelIsConjunction" ->
    MatchQ[tenGyrGeometryCalibrationHazardRel[tClaimYears, tCalibratedYears], _And],
  "absoluteRestReferenceRelIsDisjunction" ->
    MatchQ[absoluteRestReferenceRel[translationalVelocity, correctedVelocity], _Or],
  "cepheidQuasiGeometricDistanceRelIsConjunction" ->
    MatchQ[
      cepheidQuasiGeometricDistanceRel[
        distance,
        angularDiameterChange,
        radialVelocityCurve,
        pulsationCycle
      ],
      _And
    ],
  "eModePolarizationRelIsConjunction" ->
    MatchQ[eModePolarizationRel[ell, deltaTE, cEff], _And],
  "bModePolarizationRelIsConjunction" ->
    MatchQ[bModePolarizationRel[ell, deltaTB, cEff, bOverERatio], _And],
  "bOverERatioRelIsEquation" ->
    MatchQ[bOverERatioRel[bOverERatio, deltaAleph, lightSphereRadius], _Equal],
  "gutcpLightBundleConstraintRowsPresent" ->
    Length[gutcpLightBundleConstraintRows] >= 4,
  "gutcpLightBundleRelIsIntentionallyOpen" ->
    DownValues[gutcpLightBundleRel] === {},
  "readinessRows" -> readinessRows
|>;

relationSmokeReport
