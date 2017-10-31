(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
//#I "../../bin"
#r "../../src/BioFSharp.Mz/bin/Release/MathNet.Numerics.dll"
#r "../../src/BioFSharp.Mz/bin/Release/MathNet.Numerics.FSharp.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.IO.dll"
#r "../../src/BioFSharp.Mz/bin/Release/BioFSharp.Mz.dll"
//#r "../../packages/build/FSharp.Plotly/lib/net40/Fsharp.Plotly.dll"
(**
Spectrum centroidization
========================
*)

(**
This part of the documentation aims to give a brief overview of the workflow used to detect the spectral centroids of MS Spectra.
*)
open BioFSharp
open BioFSharp.Mz
open BioFSharp.IO
//open FSharp.Plotly

/// Returns the first entry of a examplary mgf File
let ms1DataTest = 
    Mgf.readMgf (__SOURCE_DIRECTORY__ + "/data/ms1Example.mgf")  
    |> List.head


#time
/// Returns a tuple of float arrays (mzData[]*intensityData[]) each containing the processed data
///

//let pickPeaksMS1 mzData intensityData = SignalDetection.Wavelet.toCentroidWithRicker2D 2 1. 0.1 0. 95. 1. mzData intensityData
//let data = SignalDetection.windowToCentroidBy pickPeaksMS1 ms1DataTest.Mass ms1DataTest.Intensity 7.5 583.3

//ms1DataTest.Mass.Length
// operation Parameters
let ms1PeakPicking mzData intensityData = 
    let parameters = SignalDetection.Wavelet.createWaveletParameters 3 1. 0.1 95. 1.
    BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D parameters mzData intensityData
 
let centroidMS1Spectrum = 
            ms1PeakPicking ms1DataTest.Mass ms1DataTest.Intensity

[
    Chart.Point( ms1DataTest.Mass, ms1DataTest.Intensity ,Name="raw data");
    Chart.Point(fst centroidMS1Spectrum,snd centroidMS1Spectrum,Name="centroided Data")    
]
|> Chart.Combine
|> Chart.withSize(500.,500.)
|> Chart.Show

#time

let ms2DataTest = 
    Mgf.readMgf (__SOURCE_DIRECTORY__ + "/data/ms2Example.mgf")  
    |> List.head


let ms2PeakPicking (mzData:float []) (intensityData: float []) = 
    if mzData.Length < 3 then 
        [||],[||]
    else
    let yThreshold = Array.min intensityData
    let paddedMz,paddedIntensity = 
        SignalDetection.Padding.paddDataWith yThreshold (Some 5) 0.05 150 50. mzData intensityData
    let parameters = SignalDetection.Wavelet.createWaveletParameters 10 yThreshold 0.1 95. 1.
    BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D parameters paddedMz paddedIntensity

let centroidMS2Spectrum = 
            ms2PeakPicking ms2DataTest.Mass ms2DataTest.Intensity
       
[
    Chart.Point(ms2DataTest.Mass, ms2DataTest.Intensity, Name="raw data");
    Chart.Point(fst centroidMS2Spectrum,snd centroidMS2Spectrum,Name="centroided Data")    
]
|> Chart.Combine
|> Chart.withSize(500.,500.)
|> Chart.Show
#time
let centroidMS2Spectrum2,corr = 
    SignalDetection.Wavelet.toCentroidWithRicker2D 10 10. 0.1 0. 90. 1. paddmzsma paddIntsma

let charts =     
    [
    for i=0 to 9 do 
        let out = ResizeArray<float>()
        corr.[i,*] |> Array.iteri (fun i x -> if i % 2 = 0 then out.Add x )
        yield out
    ] 
    |> List.map (fun x -> Chart.Point( paddmz, x))
    |> List.append [ Chart.Point( ms2DataTest.Mass, ms2DataTest.Intensity ,Name="raw data");Chart.Point(fst centroidMS2Spectrum2,snd centroidMS2Spectrum2,Name="processed data");]
    |> Chart.Combine
    |> Chart.Show
//Async.Parallel [for i = 1 to 10 do yield async {
////let centroidMS1Spectrum2 = 
//     SignalDetection.Wavelet.toCentroidWithRicker2D 10 1. 0.1 50. 95. ms1DataTest.Mass ms1DataTest.Intensity }] |> Async.RunSynchronously
(*** define-output:spectrum1 ***)

/// Creates point charts of the raw and the processed data
[
    Chart.Point( ms2DataTest.Mass, ms2DataTest.Intensity ,Name="raw data");
    Chart.Point(fst centroidMS1Spectrum2,snd centroidMS1Spectrum2,Name="raw data")
]
|> Chart.Combine
|> Chart.withSize(500.,500.)
|> Chart.Show

open MathNet.Numerics
///
let ricker2d (mzDataPadded:float []) (waveletData:float []) (centerIdx:int) (nPointsLeft:int) (nPointsRight:int) (focusedPaddedMzValue:float) (width:float) =  //(PadMz, paddedCol, nPointsLeft, nPointsRight, param1, param2, PadMz[paddedCol], waveletData)
    let inline computation cnt i =
        let vec = mzDataPadded.[i]-focusedPaddedMzValue
        let tsq = vec*vec
        let modi = 1. - tsq / (width * width)
        let gauss = exp(-1. * tsq / (2.0 * (width * width)) )
        waveletData.[cnt] <- (  2.0 / ( sqrt(3.0 * width) * (sqrt(sqrt(3.141519))) ) ) * modi * gauss
    let rec inline ricker cnt i = 
        if i = centerIdx+nPointsRight then 
                computation cnt i
        else computation cnt i
             ricker (cnt+1) (i+1)
    if nPointsLeft - nPointsRight >= waveletData.Length 
        then printfn "invalid input Parameters for ricker2d"
    ricker 0 (centerIdx-nPointsLeft)
    waveletData 

let ricker2d' (mzDataPadded:float []) (waveletData:float []) (centerIdx:int) (nPointsLeft:int) (nPointsRight:int) (focusedPaddedMzValue:float) (width:float) =  //(PadMz, paddedCol, nPointsLeft, nPointsRight, param1, param2, PadMz[paddedCol], waveletData)
    let inline computation width2 param2 param3 cnt i =
        let vec = mzDataPadded.[i]-focusedPaddedMzValue
        let tsq = vec*vec
        let modi = 1. - tsq / width2
        let gauss = exp(-1. * tsq / param2 )
        waveletData.[cnt] <- param3 * modi * gauss
    let rec inline ricker width2 param2 param3 cnt i = 
        if i = centerIdx+nPointsRight then 
                computation width2 param2 param3 cnt i
        else computation width2 param2 param3 cnt i
             ricker width2 param2 param3 (cnt+1) (i+1)
    if nPointsLeft - nPointsRight >= waveletData.Length 
        then printfn "invalid input Parameters for ricker2d"
    let width2 = (width * width)
    let param2 = (2.0 * width2)
    let param3 = 2.0 / ( sqrt(3.0 * width) * (sqrt(sqrt(3.141519))))
    ricker width2 param2 param3 0 (centerIdx-nPointsLeft)
    waveletData 
    

let data = [|1. .. 100.|] 
let mzcorr = [|1. .. 100.|] 
Chart.Point(data.[0..40], ricker2d' data mzcorr 29 20 20 30. 10.)  
|> Chart.Show

for i = 1 to 100000 do     
    ricker2d data mzcorr 29 20 20 30. 4.  

for i = 1 to 100000 do     
    ricker2d' data mzcorr 29 20 20 30. 4.  

let x = [|1. .. 20.|]
for i = 1 to 20000 do
for i in x do     
    di.TryGetValue((i,i))

let gausParamA = 
    let tmpGauss = Array.create 3 0.
    tmpGauss.[0] <- 20. 
    tmpGauss.[1] <- 100.
    tmpGauss.[2] <- 4. 
    tmpGauss |> MathNet.Numerics.LinearAlgebra.Double.DenseVector.OfArray
    
let x = [|0. .. 200.|]
let yGauss = x |> Array.map (Quantification.Fitting.Table.gaussModel.GetFunctionValue gausParamA)
Chart.Point ( x, yGauss)
|> Chart.Show

let emgParamA = 
    let tmpGauss = Array.create 4 0.
    tmpGauss.[0] <- 20. 
    tmpGauss.[1] <- 120.
    tmpGauss.[2] <- 4. 
    tmpGauss.[3] <- 8. 
    tmpGauss |> MathNet.Numerics.LinearAlgebra.Double.DenseVector.OfArray

    
let yEmg = x |> Array.map (Quantification.Fitting.Table.emgModel.GetFunctionValue emgParamA)

let finalY = Array.map2 (fun y y2 -> y + y2) yEmg yGauss

Chart.Point ( x, finalY )
|> Chart.Show

#time
let fit = Quantification.MyQuant.quantify Quantification.MyQuant.idxOfClosestLabeledPeak 11 1. 2. 120. x finalY |> fst

let yYGM = x |> Array.map (fit.Value.SelectedModel.GetFunctionValue fit.Value.EstimatedParameters)


[
Chart.Point ( x, finalY)
Chart.Point ( x, yYGM)

]
|> Chart.Combine
|> Chart.Show


(*** include-it:spectrum1 ***)

(**
If only a window of the input data shall be processed the following functions can be used.
This can be a favourable approach if only a subgroup of the data is of interest to the user 
and computation time is a limiting factor.
*)

/// Returns a tuple of float arrays (mzData[]*intensityData[]) containing only the centroids in a
/// window of a user given width centered around a user given m/z value.
let ms1CentroidsInWindow = 
     SignalDetection.windowToCentroidBy (SignalDetection.Wavelet.toSNRFilteredCentroid 0.1 50.) ms1DataTest.Mass ms1DataTest.Intensity 7.5 643.8029052

 
(*** define-output:spectrum2 ***)
/// Creates a another combined chart of the unprocessed data and the centroided data
[
    Chart.Point(ms1DataTest.Mass, ms1DataTest.Intensity,Name="raw data");
    Chart.Point(fst ms1CentroidsInWindow, snd ms1CentroidsInWindow,Name="processed data");
]
|> Chart.Combine

(*** include-it:spectrum2 ***)

/// Returns the first entry of a examplary mgf File
let ms2DataTest = 
    Mgf.readMgf (__SOURCE_DIRECTORY__ + "/data/ms2Example.mgf")  
    |> List.head

/// Returns a tuple of float arrays (mzData[]*intensityData[]) each containing the processed data
let centroidMS2Spectrum = 
    SignalDetection.Wavelet.toCentroid 0.1 ms2DataTest.Mass ms2DataTest.Intensity

//
let snrFilteredCentroidMS2Spectrum = (SignalDetection.Wavelet.toSNRFilteredCentroid 0.01 30.)   (fst centroidMS2Spectrum) (snd centroidMS2Spectrum)     
        
(*** define-output:spectrum3 ***)
/// Creates a another combined chart of the unprocessed data and the centroided MS2 data
[
    Chart.Point(ms2DataTest.Mass, ms2DataTest.Intensity,Name="raw data");
    Chart.Point(fst centroidMS2Spectrum, snd centroidMS2Spectrum,Name="centroided data");
    Chart.Point(fst snrFilteredCentroidMS2Spectrum, snd snrFilteredCentroidMS2Spectrum,Name="centroided & filtered");
]
|> Chart.Combine
(*** include-it:spectrum3 ***)
