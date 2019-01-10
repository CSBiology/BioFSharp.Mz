(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
//#I "../../bin"
//TODO: Add FILES and LINKS to other tutorials
//#I "../../bin/BioFSharp.Mz/net47/"
#I @"C:\Users\david\Source\Repos\netCoreRepos\BioFSharp.Mz\bin\BioFSharp.Mz\net47\"
#r "BioFSharp.dll"
#r "BioFSharp.Mz.dll"
#r "FSharpAux.dll"
#r "FSharp.Stats.dll"
#r "Newtonsoft.Json.dll"
#r @"C:\Users\david\Source\Repos\netCoreRepos\BioFSharp.Mz\packages\FSharp.Plotly\lib\net47\FSharp.Plotly.dll"

//open BioFSharp.Mz.SignalDetection

open BioFSharp
open BioFSharp.Mz
open ModificationInfo
open AminoAcids
open SearchDB
open BioFSharp.Digestion
open System 
open System.Runtime.Serialization.Formatters.Binary
open System.IO
open BioFSharp.Mz.SearchDB
open BioFSharp.Mz
open FSharp.Stats
open FSharp.Plotly
open FSharpAux.Array

let g1 = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|3000.;10.;0.35;|]|> vector)
let g2 = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|2500.;13.;0.35;|]|> vector)
let g3 = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|1500.;5.;0.35 ;|]|> vector)
let g4 = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|100000.;18.;0.5 ;|]|> vector)

let xData = [|0. .. 0.3 .. 80.|] |> vector

let yData = 
    let g1 = xData |> Vector.map g1 //|> vector
    let g2 = xData |> Vector.map g2 //|> vector
    let g3 = xData |> Vector.map g3 //|> vector
    let g4 = xData |> Vector.map g4 //|> vector
    
    g1 + g2 + g3 + g4
    |> Array.ofSeq
    |> Array.map (fun x -> x +  FSharp.Stats.Distributions.Continuous.Normal.Sample 0. 200.2)
    //|> Array.map (fun x -> x + 300.)





// Step 5: find rightLiftOff
let closestLiftOffIdx stepSize labeledSndDevData peakIdx   = 
    iterUntil (fun (x:Tag<Care.Extrema,(float*float)>) -> x.Meta = Care.Extrema.Negative) stepSize (peakIdx + stepSize)  labeledSndDevData 

// Step 4I: find leftLiftOffIdx
let closestLeftLiftOffIdx labeledSndDevData peakIdx =
    closestLiftOffIdx (-1) labeledSndDevData peakIdx

// Step 5: find rightLiftOff
let closestRightLiftOffIdx labeledSndDevData peakIdx = 
    closestLiftOffIdx (+1) labeledSndDevData peakIdx

 
/// Given a noisy data set, the labled negative second derivative, the index of a putative peak and the index of the peak lift of position, the function iterates
/// in the direction given by the step parameter and returns a tuple. The first value of the tuple indicates if the peak is isolated (true indicates yes) and the second value is the 
/// index index of the determined peak end. 
let tryFindPeakEnd step (xData: float []) (yData: float []) (smoothedYData: float []) (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (closestPeakIdx: int) (closestLiftOffIdx: int option) =
    printfn "clostestPeakIdx = %i" closestPeakIdx
    ///
    let signalBorderBy (step:int) =
        if Math.Sign(step) = 1 then 
            xData.Length-1
        else 
            0
    /// Inspects the sourrounding of the peak. The function walks in the direction given by the step parameter. The function accumulates all
    /// lift offs till till the next peak or the end of the signal trace is reached. Returns the last index, the number of lift offs and a bool
    /// indicating if a flanking peak is present.
    let rec loopF step (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (currentIdx: int) (kLiftOffs: int) (hasFlankingPeak:bool) = 
        if currentIdx = signalBorderBy step then 
            currentIdx, kLiftOffs, hasFlankingPeak
        else
            match kLiftOffs with 
            | x when x >=2 -> currentIdx, kLiftOffs, hasFlankingPeak
            | _ -> 
                match labeledSndDevData.[currentIdx].Meta with
                | Care.Extrema.Positive ->  
                    currentIdx, kLiftOffs, true
                | Care.Extrema.Negative ->
                    loopF step labeledSndDevData (currentIdx+step) (kLiftOffs+1) hasFlankingPeak
                |_ -> 
                    loopF step labeledSndDevData (currentIdx+step) kLiftOffs hasFlankingPeak
    match closestLiftOffIdx with 
    | Some liftOfIdx -> 
        let (_,kLiftOffs,hasFlankingPeak) = loopF step labeledSndDevData liftOfIdx 0 false 
        // Only one Liftoff and no flanking peak indicates a isolated peak.
        if kLiftOffs = 1 && hasFlankingPeak = false then
            true, 
            match iterUntili  (fun i (y:float) -> y > smoothedYData.[i-step] || y > smoothedYData.[closestPeakIdx]) step (closestPeakIdx+step) smoothedYData with 
            | None   -> signalBorderBy step
            | Some x -> x            
        // Only one Liftoff indicates a convoluted peak, iterate backwards.            
        elif kLiftOffs = 1 then
            false, 
            match  iterUntili (fun i (y:float) -> y > smoothedYData.[i-step] || y > smoothedYData.[closestPeakIdx]) step (closestPeakIdx+step) smoothedYData with
            | None   ->  (closestPeakIdx)+1
            | Some x -> x-step
        // If more than one Liftoff between two peaks is detected, the peaks are well separated
        elif kLiftOffs > 1 then
            true,  
            match iterUntili  (fun i (y:float) -> y > smoothedYData.[i-step] || y > smoothedYData.[closestPeakIdx]) step (liftOfIdx+step) smoothedYData with 
            | None   -> signalBorderBy step
            | Some x -> x        
        else
            /// No Liftoffs detected
            false,  
            match iterUntili (fun i (y:float) ->  y > smoothedYData.[i-step] || y > smoothedYData.[closestPeakIdx]) step (closestPeakIdx+step) smoothedYData with 
            | None   -> signalBorderBy step
            | Some x -> x
    | None   ->
            false, 
            match iterUntili (fun i (y:float) -> y > smoothedYData.[i-step] || y > smoothedYData.[closestPeakIdx]) step  (closestPeakIdx+step) smoothedYData with
            | None   -> signalBorderBy step 
            | Some x -> x
          
/// Given a noisy data set, the labled negative second derivative, the index of a putative peak and the index of the peak lift of position, the function iterates
/// in the positive direction returns a tuple. The first value of the tuple indicates if the peak is isolated (true indicates yes) and the second value is the 
/// index index of the determined peak end. 
let findLeftBorderOf (xData: float []) (yData: float []) smoothedYData (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (closestPeakIdx: int) (closestLiftOffIdx: int option) =
    tryFindPeakEnd (-1) xData yData smoothedYData labeledSndDevData closestPeakIdx closestLiftOffIdx

/// Given a noisy data set, the labled negative second derivative, the index of a putative peak and the index of the peak lift of position, the function iterates
/// in the positive direction returns a tuple. The first value of the tuple indicates if the peak is isolated (true indicates yes) and the second value is the 
/// index index of the determined peak end. 
let findRightBorderOf (xData: float []) (yData: float []) smoothedYData (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (closestPeakIdx: int) (closestLiftOffIdx: int option) =
    tryFindPeakEnd (1) xData yData smoothedYData labeledSndDevData closestPeakIdx closestLiftOffIdx


///
let characterizePeak (xData: float []) (yData: float []) smoothedYData (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (peakIdx: int) = 
    let apex = Some (createPeakFeature peakIdx xData.[peakIdx] yData.[peakIdx])
    let leftLiftOffIdx = closestLeftLiftOffIdx labeledSndDevData peakIdx  
    let leftLiftOff = 
        match leftLiftOffIdx with 
        | Some i -> Some (createPeakFeature i xData.[i] yData.[i])
        | None   -> None
    let convL,leftPeakEnd  = 
        let conv,leftIdx = findLeftBorderOf xData yData smoothedYData labeledSndDevData peakIdx leftLiftOffIdx
        conv, Some (createPeakFeature leftIdx  xData.[leftIdx] yData.[leftIdx])
    let rightLiftOffIdx = closestRightLiftOffIdx labeledSndDevData peakIdx  
    let rightLiftOff = 
        match rightLiftOffIdx with 
        | Some i -> Some (createPeakFeature i xData.[i] yData.[i])
        | None   -> None
    let convR,rightPeakEnd = 
        let conv,rightIdx = findRightBorderOf xData yData smoothedYData labeledSndDevData peakIdx rightLiftOffIdx
        conv, Some (createPeakFeature rightIdx xData.[rightIdx] yData.[rightIdx])
    createIdentifiedPeak 
        apex 
        leftLiftOff
        leftPeakEnd
        rightLiftOff
        rightPeakEnd
        (not convL)
        (not convR)
        xData.[leftPeakEnd.Value.Index .. rightPeakEnd.Value.Index]
        yData.[leftPeakEnd.Value.Index .. rightPeakEnd.Value.Index]


///
let filterpeaks noiseLevel (yData:float[]) (labeledDataTmp: _ []) = 
    [|
    for i = 0 to yData.Length-1 do 
        if yData.[i] > noiseLevel && labeledDataTmp.[i].Meta=SignalDetection.Care.Extrema.Positive then 
            yield i,labeledDataTmp.[i]
    |] 
                                                                                                                                                 

///    
let getPeaks snr ws xData yData = 
    ///
    let smoothedYData = Signal.Filtering.savitzky_golay ws 3 0 1 yData |> Array.ofSeq
    ///
    let negSndDev = Signal.Filtering.savitzky_golay ws 3 2 1 yData |> Array.ofSeq |> Array.map ((*) -1.)  
    ///
    let labeledDataTmp = SignalDetection.Care.labelPeaks 0. 0. (xData |> Array.ofSeq) negSndDev 
    ///
    let noiseLevel = Seq.map2 (fun x y -> abs(x-y)) smoothedYData yData |> Seq.mean |> (*) snr
    printfn "%f" noiseLevel
    /// peaks above noiselevel
    let peaks = filterpeaks noiseLevel yData labeledDataTmp // |> Array.map (fun x -> (snd x).Data)
    peaks
    //|> Array.filter (fun (i,x) -> x.Data |> fst > 14.4 && x.Data |> fst < 22. )
    |> Array.map (fun (i,_) -> characterizePeak xData yData smoothedYData labeledDataTmp i)

#time

let peaks = 
    getPeaks 2. 5 (xData |> Array.ofSeq) yData
    |> Array.map (fun x -> Chart.Line(x.XData,x.YData))
    |> Array.toList

[
Chart.Point(xData,yData)
]
@peaks
|> Chart.Combine
|> Chart.Show                                                                                                                       



//FSharp.Stats.Testing.ChiSquareTest.compute 3 [1.5;1.;0.5] [1.5;1.;0.5]

//FSharp.Stats.Testing.ChiSquareTest.compute 3 [0.6;0.3;0.1] [0.6;0.2;0.2]

//FSharp.Stats.Testing.ChiSquareTest.compute 3 [0.6;0.3;0.1] [0.6;0.2;0.2]

//let smoothedSignal = Signal.Filtering.savitzky_golay 7 3 0 3 yData 


//[
//Chart.Line(xData,yData)
//Chart.Line(xData, smoothedSignal)
//]
//|> Chart.Combine
//|> Chart.Show                                                                                                                       



//#time                                                                                                                                            
//for i = 1 to 100000 do 
//    getPeaks 0.5 7 xData yData


