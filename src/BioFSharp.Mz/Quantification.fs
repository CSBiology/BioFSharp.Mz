namespace BioFSharp.Mz

module Quantification =

    open System


    open SignalDetection
    open FSharp.Stats
    open Fitting.NonLinearRegression

    module Integration = 

        /// Returns the estimated area beneath the data using the trapezoidal rule. Requires uniform spaced datapoints. 
        let trapezEstAreaOfUniform (xData: float[]) (yData: float[]) =
            if xData.Length <> yData.Length then
                failwith "Both input collections must be of the same size"
            // lower border of area estimation
            let a = xData.[0]
            // upper border of area estimation
            let b = xData.[xData.Length-1]

            let mutable trapezSum = 0.0
    
            for i = 1 to xData.Length-2 do
                trapezSum <- trapezSum + 2. * yData.[i]
            trapezSum <- trapezSum + yData.[0]
            trapezSum <- trapezSum + yData.[xData.Length-1]
            let area = ( (b-a) / (2. * float xData.Length) ) * trapezSum
            area    
         
        /// Returns the estimated area beneath the data using the trapezoidal rule.
        let trapezEstAreaOf (xData: float[]) (yData: float[]) =
            if xData.Length <> yData.Length then
                failwith "Both input collections must be of the same size"
            let mutable trapezSum = 0.0
            for i = 0 to xData.Length-2 do
                trapezSum <- trapezSum + (xData.[i+1] - xData.[i]) * (yData.[i+1] + yData.[i])
            let area = 0.5 * trapezSum
            area    

    module GaussEstimation =
        ///
        type GaussParams = {
            Amplitude: float
            MeanX    : float
            STD      : float
            FWHM     : float
            }

        /// 
        let createGausParams amplitude meanX std fwhm = {
            Amplitude=amplitude; MeanX=meanX; STD=std; FWHM=fwhm }

        /// Returns the value of the standard deviation
        let toSTD fwhm = fwhm / ( 2. * sqrt(2. * log 2.) )

        /// Returns the full width at half maximum  
        let toFWHM std = 2. * sqrt(2. * log 2.) * std

        /// Returns the yValue of a gauss function at a given position x.
        let gaussFunc amplitude meanX std x = 
            amplitude * exp((-1.)*((((x-meanX)**2.)/(2.*std**2.))))
  
        /// Returns the yValue of a exponentially modified gauss function at a given position.
        /// The function parameter tau represents the exponential relaxation time which is the inverse of the exponential decay parameter.
        let expModGaussFunc amplitude meanX std tau x = 
            ((amplitude*std)/tau) * sqrt(System.Math.PI/2.) * exp(1./2. * ((std/tau)**2.) - ((x-meanX)/tau)) * FSharp.Stats.SpecialFunctions.Errorfunction.Erfc((1./sqrt(2.)) * ((std/tau)-((x-meanX)/std)))

        /// Estimates the Parameters of a Gaussian function
        /// Warning: This method is sensitive to noisy data. If the noise level of the input parameters is high, smoothing of 
        /// the data is strongly recommended. 
        let caruanaAlgorithm (mzData:float []) (intensityData:float []) =
            let mzData,intensityData = 
                Array.zip mzData intensityData
                |> Array.filter (fun (_,intensity) -> intensity <> 0.)
                |> Array.unzip
                |> fun (x,y) -> Vector.ofArray x, y
            if mzData.Length < 3 || intensityData.Length < 3 then None 
            else 
            let logTransIntensityData = 
                intensityData
                |> Array.map (fun x -> 
                                if x <= 0. then log 1. 
                                else log x
                             )
                |> Vector.ofArray                          
            let polCoeff = FSharp.Stats.Fitting.LinearRegression.OrdinaryLeastSquares.Polynomial.coefficient 2 mzData logTransIntensityData 
            // f(x) = a1 + a2 * x + a3 * x**2
            let a = polCoeff.[0]
            let b = polCoeff.[1]
            let c = polCoeff.[2]
            let amplitude = exp(a-((b**2.)/(4.*c)))
            let meanX = -b/(2.*c) 
            let fwhm  = sqrt(-1./(2.*c)) * sqrt(2.*log(2.))*2.
            let std   = toSTD fwhm
            Some (createGausParams amplitude meanX std fwhm)

    
        
    module MyQuant = 



                

        // TODO: Add to FSharp.Care
        let iterateTo step xData startIdx (test: 'a -> bool) =
            let rec loop  (xData: 'a []) currentIdx =
                if currentIdx <= 0 then None
                elif currentIdx >= xData.Length-1 then None
                else                                              
                    match test xData.[currentIdx] with 
                    | true -> Some currentIdx   
                    | _               -> loop xData (currentIdx+step) 
            loop xData (startIdx+step) 

        //
        type FitBothModels = 
            | True  
            | False  

        //
        type QuantificationResult = {
            FitBothModels               : FitBothModels
            EstimatedParameters         : Vector<float>
            SelectedModel               : Model
            Area                        : float
            StandardErrorOfPrediction   : float
            //If negative: MS2 recorded prior to selected peak apex, if positive: MS2 recorded after selected peak
            DeltaScanTimePeakApex       : float
            PeakApexIntensity           : float
            }

        //
        let createQuantificationResult fitBothModels estimatedParams selectedModel area standardErrorOfPrediction deltaScanTimePeakApex peakApexIntensity = {
            FitBothModels = fitBothModels;EstimatedParameters=estimatedParams; SelectedModel = selectedModel ;Area = area ; StandardErrorOfPrediction = standardErrorOfPrediction;
                DeltaScanTimePeakApex =deltaScanTimePeakApex; PeakApexIntensity=peakApexIntensity}


        ///
        let integralOfGaussian (coeffs:Vector<float>) =
            sqrt(2.) * coeffs.[0] * coeffs.[2] * sqrt(Math.PI)

        ///
        let createEMGSolverOption = createSolverOption 0.001 0.001 10000
    
        ///
        let createGaussSolverOption = createSolverOption 0.0001 0.0001 10000

        ///
        let findRightFittingIdx (xAndYData: (float*float) []) (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (closestPeakIdx: int) (closestRightLiftOffIdx: int option) =
            let rec loopF (labeledSndDevData: Tag<Care.Extrema,(float*float)> []) (currentIdx: int) (kLiftOffs: int) (hasRightPeak:bool) = 
                if currentIdx = labeledSndDevData.Length-1 then 
                    currentIdx, kLiftOffs, hasRightPeak
                else
                    match labeledSndDevData.[currentIdx] with
                    | x when x.Meta = Care.Extrema.Positive ->  
                        currentIdx, kLiftOffs, true
                    | x when x.Meta = Care.Extrema.Negative ->
                        loopF labeledSndDevData (currentIdx+1) (kLiftOffs+1) hasRightPeak
                    |_ -> 
                        loopF labeledSndDevData (currentIdx+1) kLiftOffs hasRightPeak
            match closestRightLiftOffIdx with 
            | Some x -> 
                let (currentIdx, kLiftOffs, hasRightPeak) = loopF labeledSndDevData closestRightLiftOffIdx.Value 0 false 
                    // only one Liftoff and no flanking peak indicates a isolated peak and both models can be tested. 
                if kLiftOffs = 1 && hasRightPeak = false then
                    FitBothModels.True, 
                    match iterateTo (+1) xAndYData (closestRightLiftOffIdx.Value) (fun (x:float*float) -> snd x <= snd xAndYData.[closestRightLiftOffIdx.Value] || snd x >= snd xAndYData.[closestRightLiftOffIdx.Value]) with 
                    | None -> xAndYData.Length-1
                    | Some x -> x            
                // only one Liftoff indicates a convoluted peak, use only Gaussian model            
                elif kLiftOffs = 1 then
                    FitBothModels.False, 
                        match iterateTo (-1) xAndYData (closestRightLiftOffIdx.Value) (fun (x:float*float) -> snd x > snd xAndYData.[closestRightLiftOffIdx.Value]) with
                        | None ->  (closestPeakIdx)+1
                        | Some x -> x
                // if more than one Liftoff between two peaks is detected, the peaks are well separated and both Models can be tested
                elif kLiftOffs > 1 then 
                    FitBothModels.True,  
                    match iterateTo (+1) xAndYData (closestRightLiftOffIdx.Value) (fun (x:float*float) -> snd x <= snd xAndYData.[closestRightLiftOffIdx.Value] || snd x > snd xAndYData.[closestRightLiftOffIdx.Value]) with 
                    | None -> xAndYData.Length-1
                    | Some x -> x        
                else
                    FitBothModels.False,  
                    match iterateTo (+1) xAndYData (closestPeakIdx) (fun (x:float*float) ->  snd x < 0.5 * snd xAndYData.[closestPeakIdx]) with 
                    | None -> xAndYData.Length-1
                    | Some x -> x
            | None   -> 
                FitBothModels.False, 
                    match iterateTo (+1) xAndYData (closestPeakIdx) (fun (x:float*float) -> snd x < 0.5 * snd xAndYData.[closestPeakIdx]) with
                    | None   -> xAndYData.Length-1 
                    | Some x -> x
            
 
   
        /// 
        let quantify (peakByF:Tag<Care.Extrema,(float*float)> [] -> Care.Extrema -> 'a -> ((int * Tag<Care.Extrema,(float*float)>) option)) windowSizeSGfilter negYThreshold posYThreshold (scanTime: float) (xData :float []) (yData: float [])= 

            printfn "start quantification"
            if xData.Length < 6 || yData.Length < 6 then None
            else
            // Step 0: zip xData and yData
            let xAndYData = 
                Array.zip xData yData
            // Step 1: Calculate negative snd derivative of the intensity Data
            let negSndDev = 
                Signal.Filtering.savitzky_golay windowSizeSGfilter 3 2 3 yData 
                |> Array.ofSeq
                |> Array.map (fun x -> x * -1.)    
            // Step 2: label data points to be local Minima or maxima
            let labeledSndDevData = 
                let labeledDataTmp = SignalDetection.Care.labelPeaks negYThreshold posYThreshold xData negSndDev
                let maxPeakIntensity = 
                    labeledDataTmp 
                    |> Array.maxBy (fun x -> x.Meta = Care.Extrema.Positive)
                labeledDataTmp 
                |> Array.map (fun x -> if x.Meta = Care.Extrema.Positive then
                                         match snd x.Data with
                                         | validIntensity when validIntensity > 0.05 * (snd maxPeakIntensity.Data) -> x                             
                                         | _                                                                       -> 
                                            {Meta=Care.Extrema.None; Data= x.Data}
                                       else x
                             )
            // Step 3: find closest Peak to MS2 scantime
            let closestPeakIdx = 
                peakByF labeledSndDevData SignalDetection.Care.Extrema.Positive scanTime
            if closestPeakIdx.IsNone then None
            else
            // Step 4I: find leftLiftOffIdx
            let closestLeftLiftOffIdx =
                iterateTo (-1) labeledSndDevData (fst closestPeakIdx.Value) (fun (x:Tag<Care.Extrema,(float*float)>) -> x.Meta = Care.Extrema.Negative)
            // Step4II: find leftFittingStartIdx
            let leftfittingBorderIdx = 
                match closestLeftLiftOffIdx with 
                | None   -> 0 
                | Some x -> x+1
            // Step 5: find rightLiftOff
            let closestRightLiftOffIdx = 
                iterateTo (+1) labeledSndDevData (fst closestPeakIdx.Value) (fun (x:Tag<Care.Extrema,(float*float)>) -> x.Meta = Care.Extrema.Negative)
            // Step 6: check if another peak is present to the right side of the chosen peak, count Liftoffs points, determine Model selection     
            let (modelStatus,rightFittingBorderIdx) =
                findRightFittingIdx xAndYData labeledSndDevData (fst closestPeakIdx.Value) closestRightLiftOffIdx
            //Step 7: Create sub-array of xData and yData that is considered for subsequent fitting procedures    
            let xDataForFit = xData.[leftfittingBorderIdx.. rightFittingBorderIdx] 
            let yDataForFit = yData.[leftfittingBorderIdx.. rightFittingBorderIdx]
            //Step 8: Use Caruanas algorithm to estimate parameters height, position and Full width at half maximum of 
            //        the selected peak
            let gausParamEstCaruana = 
                GaussEstimation.caruanaAlgorithm xDataForFit yDataForFit
            if gausParamEstCaruana.IsNone then None
            else 
            let gausParamEstCaruana = gausParamEstCaruana.Value
            //Step 9: Case A: if FitBithModels = True, the peak ending can be used to estimate a possible tailing    if FitBothModels = False then the first peak of a convoluted peak pair was chosen and is subsequently used to estimate the area
            //        Case B: if FitBothModels = False then the first peak of a convoluted peak pair was chosen and is subsequently used to estimate the area
            // Case A:
            if modelStatus = FitBothModels.True then
                let modelFunction = 
                    ///
                    let gausParamA = 
                        let tmpGauss = Array.create 3 0.
                        tmpGauss.[0] <- gausParamEstCaruana.Amplitude 
                        tmpGauss.[1] <- (gausParamEstCaruana.MeanX ) 
                        tmpGauss.[2] <- gausParamEstCaruana.STD 
                        tmpGauss
                    ///
                    let gaussPrediction =
                        let gaussSolOptions = createGaussSolverOption gausParamA
                        try
                            let gaussParamA = 
                                    //    let lambdaInitial = 0.001
        //    let lambdaFactor  = 10.0
                                LevenbergMarquardt.estimatedParams Table.gaussModel gaussSolOptions 0.001 10.0 xDataForFit yDataForFit
                            //FShapFitting.levenbergMarquardtSolver Fitting.Table.gaussModel gaussSolOptions xDataForFit yDataForFit paramConti                                   
                            let gaussYPredicted = Array.map (fun xValue -> Table.gaussModel.GetFunctionValue gaussParamA xValue) xDataForFit
                            Some (gaussParamA, gaussYPredicted)
                        with 
                        | _ as ex  -> 
                            None
          
                    ///
                    let exponentialDecayEst = 
                        let startTime = xData.[closestRightLiftOffIdx.Value]
                        let startIntensity = yData.[closestRightLiftOffIdx.Value]
                        let idxHalfIntensity = iterateTo (+1) yData closestRightLiftOffIdx.Value (fun x -> x < 0.5*startIntensity)    
                        let endTime = 
                            match idxHalfIntensity with
                            | Some x -> xData.[x]
                            | None   -> xData.[xData.Length-1]
                            
                        let decayEst = (endTime-startTime) / (log 2.)
                        decayEst
                    
                    ///
                    let gausParamEMG = 
                        let tmpEmg = Array.create 4 0.
                        tmpEmg.[0] <- gausParamEstCaruana.Amplitude 
                        tmpEmg.[1] <- (gausParamEstCaruana.MeanX )
                        tmpEmg.[2] <- gausParamEstCaruana.STD 
                        tmpEmg.[3] <- exponentialDecayEst
                        tmpEmg
 
                    ///
                    let emgPrediction =
                        [|exponentialDecayEst*0.75 ..  exponentialDecayEst|]
                        |> Array.map (fun x -> 
                                        gausParamEMG.[3] <- x
                                        let emgSolOptions = createEMGSolverOption gausParamEMG
                                        try 
                                            let emgParamA = 
                                                LevenbergMarquardt.estimatedParams Table.emgModel emgSolOptions 0.001 10.0 xDataForFit yDataForFit
                                                //Fitting.levenbergMarquardtSolver Fitting.Table.emgModel emgSolOptions xDataForFit yDataForFit paramConti
                                            let emgYPredicted = Array.map (fun xValue -> Table.emgModel.GetFunctionValue emgParamA xValue) xDataForFit
                                            Some (emgParamA, emgYPredicted)
                                        with 
                                        | _ as ex  -> 
                                                None
                                    )                
                        |> Array.filter (fun x -> x.IsSome)
                        |> fun possibleEMGFits -> 
                            match possibleEMGFits with
                            | pEmgF when Array.isEmpty pEmgF -> None
                            | _                      ->
                                possibleEMGFits  
                                |> Array.map (fun x ->
                                                let sEoE = standardErrorOfPrediction 4. (snd x.Value) yDataForFit  
                                                x, sEoE
                                            )
                                |> Array.minBy (fun x -> snd x)
                                |> fun x -> fst x
            
                    // If both models are fittable, choose the model with smaller standard error of the estimate
                    if gaussPrediction.IsSome && emgPrediction.IsSome then
                        ///
                        let yGauss = snd gaussPrediction.Value
                        let sEoE_Gauss = standardErrorOfPrediction 3. yGauss yDataForFit
                        ///
                        let yEMG   = snd emgPrediction.Value
                        let sEoE_EMG   = standardErrorOfPrediction 4. yEMG yDataForFit 
                        if sEoE_Gauss > -1. && sEoE_EMG > -1. && sEoE_Gauss > sEoE_EMG then
                            Some (Table.emgModel, emgPrediction.Value, sEoE_EMG)
                        else        
                            Some (Table.gaussModel, gaussPrediction.Value, sEoE_Gauss)
                    elif emgPrediction.IsSome then
                        let yEMG   = snd emgPrediction.Value
                        let sEoE_EMG   = standardErrorOfPrediction 4. yEMG yDataForFit      
                        Some (Table.emgModel , emgPrediction.Value, sEoE_EMG)
                    elif gaussPrediction.IsSome then
                        let yGauss = snd gaussPrediction.Value
                        let sEoE_Gauss = standardErrorOfPrediction 3. yGauss yDataForFit
                        Some (Table.gaussModel, gaussPrediction.Value, sEoE_Gauss)
                    else
                        None
                // compute area beneath curve
                match modelFunction with
                | None -> 
                    None       
                | Some (modelF,(paramV, yData), sEoE) when paramV.Length = 3 ->
                    try 
                        let area = 
                            //MathNet.Numerics.Integrate.OnClosedInterval((fun x -> (modelF.GetFunctionValue paramV x)),xData.[0],xData.[xData.Length-1])
                            integralOfGaussian paramV
                        let deltaScanTimePeakApex = (scanTime - gausParamEstCaruana.MeanX)
                        //
                        Some (createQuantificationResult FitBothModels.True paramV modelF area sEoE deltaScanTimePeakApex (modelF.GetFunctionValue paramV gausParamEstCaruana.MeanX))
                    with 
                    | _ as ex ->    None
                | Some (modelF,(paramV, yData), sEoE) when paramV.Length = 4 ->
                    try 
                        let f = modelF.GetFunctionValue paramV
                        let integrationXData =
                            [|xData.[0] .. 0.02 .. xData.[xData.Length-1]|]
                        let integrationYData = 
                            integrationXData |> Array.map f
                        let area = 
                            Integration.trapezEstAreaOf integrationXData integrationYData
                            //MathNet.Numerics.Integrate.OnClosedInterval((fun x -> (modelF.GetFunctionValue paramV x)),xData.[0],xData.[xData.Length-1])
                        let deltaScanTimePeakApex = (scanTime - gausParamEstCaruana.MeanX)
                        //
                        Some (createQuantificationResult FitBothModels.True paramV modelF area sEoE deltaScanTimePeakApex (modelF.GetFunctionValue paramV gausParamEstCaruana.MeanX))
                    with 
                    | _ as ex ->    None
            // Case B:
            else
                let modelFunction = 
                    ///
                    let gausParamA = 
                        let tmpGauss = Array.create 3 0.
                        tmpGauss.[0] <- gausParamEstCaruana.Amplitude 
                        tmpGauss.[1] <- (gausParamEstCaruana.MeanX ) 
                        tmpGauss.[2] <- gausParamEstCaruana.STD 
                        tmpGauss
                    ///
                    let gaussPrediction =
                        let gaussSolOptions = createGaussSolverOption gausParamA
                        try
                            let gaussParamA = LevenbergMarquardt.estimatedParams Table.gaussModel gaussSolOptions 0.001 10.0 xDataForFit yDataForFit                                                    
                            let gaussYPredicted = Array.map (fun xValue -> Table.gaussModel.GetFunctionValue gaussParamA xValue) xDataForFit
                            Some (gaussParamA, gaussYPredicted)
                        with 
                        | :? System.ArgumentException as ex  -> 
                            None
                       
                    if gaussPrediction.IsSome then
                        let yGauss = snd gaussPrediction.Value
                        let sEoE_Gauss = standardErrorOfPrediction 3. yGauss yDataForFit
                        Some (Table.gaussModel, gaussPrediction.Value, sEoE_Gauss)
                    else
                        None
                // compute area beneath curve
                match modelFunction with
                | None -> 
                    None      
                | Some (modelF,(paramV, yData), sEoE) ->
                    try                              
                        //
                        let area = 
                            //MathNet.Numerics.Integrate.OnClosedInterval((fun x -> (modelF.GetFunctionValue paramV x)),xData.[0],xData.[xData.Length-1]) 
                            integralOfGaussian paramV
                        //
                        let deltaScanTimePeakApex = scanTime - gausParamEstCaruana.MeanX
                        Some (createQuantificationResult FitBothModels.False paramV modelF area sEoE deltaScanTimePeakApex (modelF.GetFunctionValue paramV gausParamEstCaruana.MeanX))
                    with 
                    | _ as ex -> 
                        None

                    