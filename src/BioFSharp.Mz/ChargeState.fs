﻿namespace BioFSharp.Mz


module ChargeState =
    
    open System
    open BioFSharp
    open FSharpAux

    open Peaks
    open PeakArray
    open PeakList
    open SignalDetection

    
    type ChargeDetermParams = {
        ExpectedMinimalCharge : int ///TODO: learn from Data
        ExpectedMaximumCharge : int ///TODO: learn from Data
        Width: float
        /// RelativeToStartPeak
        MinIntensity: float
        /// RelativeToPriorPeak
        DeltaMinIntensity: float
        NrOfRndSpectra        : int
        }

    let createChargeDetermParams expectedMinimalCharge expectedMaximumCharge width minIntensity deltaMinIntensity nrOfRndSpectra = {
        ExpectedMinimalCharge=expectedMinimalCharge; ExpectedMaximumCharge=expectedMaximumCharge; Width=width; MinIntensity=minIntensity; DeltaMinIntensity=deltaMinIntensity; NrOfRndSpectra=nrOfRndSpectra }
    
    type PutativeIsotopeCluster = {
        Peaks           : PeakList<Peak>   
        SourceSetLength : int 
        SubSetLength    : int
        }

    let createPutativeIsotopeCluster peaks sourceSetLength subSetLength = {
        Peaks=peaks; SourceSetLength=sourceSetLength; SubSetLength=subSetLength
        }

    type AssignedCharge = { 
        PrecMZ                      : float
        Charge                      : int
        PutMass                     : float
        MZChargeDev                 : float
        Score                       : float
        DistanceRealTheoPeakSpacing : float list
        SubSetLength                : int 
        StartPeakIntensity          : float
        Peaks        : Set<Peak>
        } 

    let createAssignedCharge precMZ charge putMass mZChargeDev score distanceRealTheoPeakSpacing subSetLength startPeakIntensity peaks= {
        PrecMZ=precMZ; Charge=charge; PutMass=putMass; MZChargeDev=mZChargeDev; Score=score;DistanceRealTheoPeakSpacing=distanceRealTheoPeakSpacing; 
            SubSetLength=subSetLength; StartPeakIntensity=startPeakIntensity; Peaks=peaks }
   
    type TestedItem<'a> = {
        TestedObject: 'a
        PValue: float
        }     
    
    let createTestedItem testedObject pValue = {
        TestedObject=testedObject; PValue=pValue }

    /// Returns a Collection of MZIntensityPeaks, The Collection starts with the first Element on the right side of the startIdx. 
    /// and ends either with the last element of the mzIntensityArray or when the MzDistance to the highest Peak exceeds 
    /// the given windowwidth.   
    let getRelPeakPosInWindowBy (mzData: float []) (intensityData: float []) width minIntensity deltaMinIntensity startIdx =
        let mzDataLength = mzData.Length
        if mzDataLength > 0 then
            let startPkMZ = mzData.[startIdx]
            let startPkInt = intensityData.[startIdx]
            let hasValidIntensity count =
                intensityData.[count] > minIntensity * startPkInt &&
                intensityData.[count] > deltaMinIntensity * intensityData.[count-1]   
            let rec loop count accmz (mzData: float []) (intensityData: float []) =            
                if count = mzDataLength then 
                     accmz
                elif mzData.[count] - startPkMZ > width then 
                     accmz
                elif count = mzDataLength-1 && hasValidIntensity count then 
                     loop (count+1) ((createPeak (mzData.[count]-startPkMZ) (intensityData.[count] / startPkInt) )::accmz) mzData intensityData
                elif hasValidIntensity count then 
                     loop (count+1) ((createPeak (mzData.[count]-startPkMZ) (intensityData.[count] / startPkInt) )::accmz) mzData intensityData
                else loop (count+1) accmz mzData intensityData
            loop (startIdx+1) [] mzData intensityData
            |> fun sourceSet -> startPkInt , createPutativeIsotopeCluster sourceSet (sourceSet.Length+1) (sourceSet.Length+1) //the length must be raised by 1 because the first element (0.,1.) is left out from the collection 
        else 0., createPutativeIsotopeCluster [] 0 0
    
    /// Creates the PowerSet of a given Collection of MZIntensityPeaks. Adds a StartPeak with the relative Position 0 and the 
    /// relative Intensity 1 to each subSet. 
    let powerSetOf (putCluster: PutativeIsotopeCluster) = 
        let rec createPowerset (inputL: PeakList<_>) = 
            match inputL with
            | [] -> [[createPeak 0. 1.]] 
            | x::xs -> List.collect (fun subSet ->
                                                    [subSet; x::subSet]
                                    ) (createPowerset xs)
        let superSet = createPowerset (putCluster.Peaks) 
        superSet
        |> List.map (fun subSet -> createPutativeIsotopeCluster subSet putCluster.SourceSetLength subSet.Length)
   
    /// Calculates the mzDistances of a List of MzIntensityPeaks
    let mzDistancesOf  (mzIntensityPeaks: PeakList<_>) = //TODO: calc Slope; times Slope changed
        let rec innerF acc (mzIntensityPeaks: PeakList<_>) =
            match mzIntensityPeaks with 
            | []      -> acc
            | h::[]   -> acc 
            | h::tail -> innerF (h.Mz - tail.Head.Mz::acc) tail
        innerF [] mzIntensityPeaks

    /// Returns a charge state based on a given meanOfInterPeakDistances. The list of possible charge states is defined by
    /// the user.  
    let getChargeBy (chargeDeterminationParams: ChargeDetermParams) (meanOfInterPeakDistances: float) = 
        let chargeL = [chargeDeterminationParams.ExpectedMinimalCharge.. chargeDeterminationParams.ExpectedMaximumCharge]
        let rec innerF (acc:int*float) (chargeL: int list) (meanOfPeakDistances:float) = 
                match chargeL with 
                | [] -> fst acc
                | h::tail -> let theoreticalInterIsotopeDistance = 1./ float h
                             let difMeanPeakDistanceToTheoDistance = meanOfPeakDistances-theoreticalInterIsotopeDistance |> abs
                             if theoreticalInterIsotopeDistance <= meanOfPeakDistances && difMeanPeakDistanceToTheoDistance < snd acc then 
                                  h  
                             elif difMeanPeakDistanceToTheoDistance > snd acc then 
                                  (fst acc)
                             else innerF (h, difMeanPeakDistanceToTheoDistance ) tail meanOfPeakDistances 
        innerF (0, 100.) chargeL meanOfInterPeakDistances 

    /// Returns the MZChargeDeviation based on the theoreticalInterIsotopeDistance as a Measure for central tendency
    /// at a given ChargeState. 
    let mzChargeDeviationBy (interPeakDistances:seq<float>) (theoreticalInterIsotopeDistance)  =
        use e = interPeakDistances.GetEnumerator()
        let rec loop n (acc:float) =
            match e.MoveNext() with
            | true  -> loop (n + 1) (acc + ((square (e.Current - theoreticalInterIsotopeDistance) )/theoreticalInterIsotopeDistance))
            | false -> if (n > 0) then sqrt (acc / (float n)) else nan            
        loop 0 0.0 

    /// Returns a empirically determined Score to achieve a Ranking of SubSets. The optimal weighting of the parameters was
    /// determined via Linear Discriminant Analysis.  
    let getScore (elementsInSubSet:int) (elementsInSourceSet:int) deviationPeakDistanceToTheoDistance = 
        let ldaVec  = [|1.0; -0.30|]
        ldaVec.[0] * deviationPeakDistanceToTheoDistance + ldaVec.[1] * (float elementsInSubSet/ float elementsInSourceSet)  

    /// Returns a random integer between a lower and a upper Value
    let rndIntBetween (rnd:System.Random) lowerBorder upperBorder  = rnd.Next(lowerBorder, upperBorder)

    /// Returns a possible InterPeakDistance by a given charge state. The retrieved distance follows a normaldistribution
    /// centered around the theoretical interPeakDistance. The standardDeviation is dependent on the used mass spectrometer
    let interPeakDistanceBy (stdDev: float)  (assignedCharge: float) = 
        let calcRndPosition assignedCharge stdDev =
            1./assignedCharge + (FSharp.Stats.Distributions.Continuous.normal 0. stdDev).Sample()
        calcRndPosition assignedCharge stdDev

    /// Creates a random MzIntensityEntityCollection  
    let rndMzIntensityEntityCollectionBy (rnd:System.Random) (stdDev: float) (maxCharge: int) maxDistance minChainLength =
        let rec innerF count (accFloat:float) (accList:PeakList<_>) minChainLength  =
            let rndCharge = 
                rndIntBetween rnd 1 maxCharge
                |> float 
            let nextRnd = (interPeakDistanceBy stdDev rndCharge)
            if (accFloat+nextRnd) > maxDistance then 
                 innerF 2 0. [createPeak 0. 1.] (minChainLength)  
            elif count = minChainLength then 
                 createPutativeIsotopeCluster ((createPeak (accFloat+nextRnd) 0.)::accList) minChainLength minChainLength    
            else innerF (count+1) (accFloat+nextRnd) ((createPeak (accFloat+nextRnd) 0.)::accList) minChainLength     
        innerF 2 0. [createPeak 0. 1.] (minChainLength)

        
    /// Creates a user defined amount of random spectra of defined length. Returns the mzChargeDeviation of each simulated Spectrum
    let mzDevOfRndSpec (rnd:System.Random) chargeStateDetermParams peakPosStdDev (chainLength,charge)  = 
        [1.. chargeStateDetermParams.NrOfRndSpectra] 
        |> List.map (fun i -> rndMzIntensityEntityCollectionBy (rnd:System.Random) peakPosStdDev chargeStateDetermParams.ExpectedMaximumCharge chargeStateDetermParams.Width chainLength)           
        |> List.map (fun subSet -> 
                        let interPeakDistances = mzDistancesOf subSet.Peaks                    
                        let mzChargeDeviation = mzChargeDeviationBy interPeakDistances (1./charge)  
                        mzChargeDeviation
                        )
        |> List.sort
        |> Array.ofList
    
    /// Returns Function to generate random spectra and to calculate their mzChargeDeviations.   
    let initMzDevOfRndSpec (rnd:System.Random) (chargeStateDetermParams: ChargeDetermParams) peakPosStdDev =
        Memoization.memoize (mzDevOfRndSpec rnd (chargeStateDetermParams: ChargeDetermParams) peakPosStdDev) ///stdv 0.01516580549
        
    /// Returns the empirically determined PValue. The PValue is the quotient of simulated mzChargeDeviations lower than the mzChargeDeviation
    /// observed divided by their total number
    let empiricalPValueOfSim initGenerateMzSpecDevWithMemF (nrOfPeaksInSubSet,charge) score  = //TODO nrOfPeaks,charge score in parameter
        let generateMzSpecDev = initGenerateMzSpecDevWithMemF (nrOfPeaksInSubSet,charge)
        let numerator =  (generateMzSpecDev |> Array.tryFindIndex (fun x -> x > score)) 
        match numerator with
        | Some x -> (float x) / float generateMzSpecDev.Length
        | None -> 1.

    /// Returns the empirically determined PValue. The PValue is the quotient of simulated mzChargeDeviations lower than the mzChargeDeviation
    /// observed divided by their total number
    let empiricalRightPValueOf distribution score  = //TODO nrOfPeaks,charge score in parameter
        let numerator =  (distribution |> Array.tryFindIndex (fun x -> x > score)) 
        match numerator with
        | Some x -> (float x) / float distribution.Length
        | None -> 1.

    /// Returns list of putative precursorChargeStates along with Properties used for evaluation.
    let putativePrecursorChargeStatesBy (chargeDeterminationParams: ChargeDetermParams) (mzData: float []) (intensityData: float []) (precursorMZ:float) =
        if mzData |> Array.isEmpty || intensityData |> Array.isEmpty then []
        else
        let (startPeakIntensity,originSet) = getRelPeakPosInWindowBy (mzData: float []) (intensityData: float [])  chargeDeterminationParams.Width chargeDeterminationParams.MinIntensity chargeDeterminationParams.DeltaMinIntensity (Care.idxOfClosestPeakBy  (mzData: float []) (intensityData: float [])  precursorMZ)
        originSet
        |> powerSetOf 
        |> List.filter (fun subSet -> subSet.SubSetLength > 1)
        |> List.map (fun subSet -> 
                        let peakPos = 
                            subSet.Peaks
                            |> Set.ofList
                        let interPeakDistances = mzDistancesOf subSet.Peaks 
                        let meanInterPeakDistances = FSharp.Stats.List.mean interPeakDistances
                        let assignedCharge = getChargeBy chargeDeterminationParams meanInterPeakDistances
                        let theoInterPeakDistances = 1. / float assignedCharge
                        let distanceRealTheoPeakSpacing = 
                            interPeakDistances
                            |> List.map (fun distance -> distance - theoInterPeakDistances ) 
                        let mzChargeDeviation = mzChargeDeviationBy interPeakDistances theoInterPeakDistances
                        let score = getScore subSet.SubSetLength subSet.SourceSetLength mzChargeDeviation
                        let putMass = BioFSharp.Mass.ofMZ precursorMZ (float assignedCharge)
                        createAssignedCharge precursorMZ assignedCharge putMass mzChargeDeviation score distanceRealTheoPeakSpacing subSet.SubSetLength startPeakIntensity peakPos  
                     )
        |> List.sortBy (fun assignedCharge ->  assignedCharge.Score)
        |> List.distinctBy (fun assignedCharge ->  assignedCharge.Charge)
                                     
    /// Returns the StandardDeviation of the PeakDistances
    let peakPosStdDevBy (putativeChargeStates: AssignedCharge list) = 
        putativeChargeStates
        |> List.map (fun putativeChS -> putativeChS.DistanceRealTheoPeakSpacing )
        |> List.concat 
        |> FSharp.Stats.Seq.stDev

    /// Returns a List of tested AssignedCharges. This Function eliminates all
    let removeSubSetsOfBestHit (assignedCharges: TestedItem<AssignedCharge> list) =
        let rec loop bestSet acc (assignedCharges: TestedItem<AssignedCharge> list) =
            match assignedCharges with 
            | [] -> (bestSet::acc) |> List.sortBy (fun x -> x.TestedObject.Score) 
            | h::tail -> match Set.isSuperset bestSet.TestedObject.Peaks h.TestedObject.Peaks with
                         | false  -> loop bestSet (h::acc) tail 
                         | true -> loop bestSet acc tail
        if assignedCharges = [] then 
            [] 
        else
            let bestSet = assignedCharges.Head
            loop bestSet [] assignedCharges

    ///
    let normalizePeaksByIntensitySum (peaks: Set<BioFSharp.Mz.Peak>) =
        let tmp = peaks |> Set.toArray
        let sumOfIntensities = 
            tmp |> Array.sumBy (fun x -> x.Intensity)        
        tmp
        |> Array.map (fun x -> x.Intensity / sumOfIntensities)

    ///
    let n14MassToLambda mass =  mass * 0.0005903277935 - 0.03542317929

    ///
    let n15MassToLambda mass =  mass * 0.0005306674223 - 0.03184171518
    
    ///
    let poissonProb lambda xValue =
        ( ( lambda ** xValue) / FSharp.Stats.SpecialFunctions.Factorial.factorial (int xValue) ) * exp(-lambda) 

    ///
    let poissonEstofMassTrunc (massToLambda: float -> float) (limit:int) mass = 
        let tmp = 
            [|for i=0. to float limit-1. do 
                 yield poissonProb (massToLambda mass) i
            |]
        let SumOfProb = 
            tmp |> Array.sum
        tmp 
        |> Array.map (fun x -> x / SumOfProb)
         
    ///
    let kullbackLeiblerDivergenceOf (qTo:float []) (p:float []) = 
        if qTo.Length <> p.Length then None
        else
        let divergence =
            Array.fold2 (fun acc qY pY -> 
                            qY
                            |> (/) pY
                            |> log
                            |> (*) pY
                            |> (+) acc
                        ) 0. qTo p 
        Some divergence
