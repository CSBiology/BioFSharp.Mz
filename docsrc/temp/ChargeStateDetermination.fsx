(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
//#I "../../bin"
#r "../../src/BioFSharp.Mz/bin/Release/Fsharp.Care.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.Mz.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.IO.dll"
#r "../../src/BioFSharp.Mz/bin/Release/MathNet.Numerics.dll"
#r "../../src/BioFSharp.Mz/bin/Release/MathNet.Numerics.FSharp.dll"
#r "../../packages/build/FSharp.Plotly/lib/net40/Fsharp.Plotly.dll"
(**
Charge state determination
==========================

The charge state determination algorithm implemented in ChargeState.fs analyzes the spacing in between 
single peaks of a isotope Cluster. 
*)

(**
Returns the first entry of a examplary mgf File
*)

open System
open BioFSharp
open BioFSharp.IO
open BioFSharp.Mz
open FSharp.Care
open FSharp.Plotly
open MathNet.Numerics


///Returns the first entry of a examplary mgf File
let ms1DataTest = 
    Mgf.readMgf (__SOURCE_DIRECTORY__ + "/data/ms1Example.mgf")  
    |> List.head

(**
TestParameters used for chargestate determination
*)

///TestParameters used for chargestate determination
let chargeDetParamTest  = 
    ChargeState.createChargeDetermParams
        // minimum expected chargestate
        1
        // maximum expected chargestate
        8
        // width of the peakwindow used for chargestate determination
        1.1
        // minimal relative peak intensity in respect to the initial peak in the peakwindow
        0.15
        // minimal relative peak intensity in respect to the previous peak in the peakwindow
        0.3
        // number of theoretical spectra computed used to determine the pValue
        10000

(**
PrecursorMZ that was picked in the MS1 full scan to generate a MS2
*)

/// PrecursorMZ that was picked in the MS1 full scan to generate a MS2
let ms2PrecursorMZ =  643.8029052

(**
Instance of a Mersenne Twister, a random number generator used to compute the position of a peak in a 
random spectrum
*)

/// Instance of a Mersenne Twister, a random number generator used to compute the position of a peak in a 
/// random spectrum
let rnd = new System.Random()

(**
Returns a function that generates random spectra of a given length and calculates the mzDeviation assuming
a given chargestate
*)

/// Returns a function that generates random spectra of a given length and calculates the mzDeviation assuming
/// a given chargestate
let initGen = ChargeState.initMzDevOfRndSpec rnd chargeDetParamTest 0.015//peakPosStdDev

(**
Returns a tuple of float arrays (mzData[]*intensityData[]) containing only the centroids in a
window of a user given width centered around a user given m/z value.
*)

/// Returns a tuple of float arrays (mzData[]*intensityData[]) containing only the centroids in a
/// window of a user given width centered around a user given m/z value.
let centroidsInWindow = 
    SignalDetection.windowToCentroidBy (SignalDetection.Wavelet.toCentroidWithRicker2D 10 1. 0.1 50. 95.) ms1DataTest.Mass ms1DataTest.Intensity 7.5 ms2PrecursorMZ  

(**
Returns a list of assigned chargestates sorted by their score.
*)

/// Returns a list of assigned chargestates sorted by their score.
let putativeCharges = 
    ChargeState.putativePrecursorChargeStatesBy chargeDetParamTest (fst centroidsInWindow) (snd centroidsInWindow) ms2PrecursorMZ

(**
Returns a list of assigned chargestates after hypothesis testing. This list is filtered for items that meet the 
significance level
*)
/// Returns a list of assigned chargestates after hypothesis testing. This list is filtered for items that meet the 
/// significance level
let testedItems = 
    putativeCharges 
    |> List.map (fun assCh -> ChargeState.createTestedItem assCh (ChargeState.empiricalPValueOfSim initGen (assCh.SubSetLength ,float assCh.Charge) assCh.MZChargeDev ))
    |> ChargeState.removeSubSetsOfBestHit
    |> List.filter (fun item -> item.PValue < 0.05)
    |> List.map (fun testedI -> testedI.TestedObject)
    |> List.map (fun testedI -> 
                    let normPeaks = 
                        ChargeState.normalizePeaksByIntensitySum testedI.Peaks 
                    let poissEst = 
                            ChargeState.poissonEstofMassTrunc ChargeState.n15MassToLambda normPeaks.Length testedI.PutMass
                    ChargeState.kullbackLeiblerDivergenceOf normPeaks poissEst
                )
                        

