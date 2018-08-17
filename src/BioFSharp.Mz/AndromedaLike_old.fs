namespace BioFSharp.Mz

open System
open BioFSharp
open BioFSharp.Mz.Peaks
open SearchDB
open Fragmentation
open TheoreticalSpectra
open SearchEngineResult

open MathNet.Numerics

module AndromedaLike2 =

    1

    //type MatchingScore = {
    //    Score   : float 
    //    N       : int 
    //    K       : int
    //    Q       : int
    //    //Matches : PeakAnnotation list 
    //    }
    
    //let createMatchingScore score n k q matches = {
    // Score = score; N = n ; K = k; Q = q ;//Matches = matches
    //}

    /////
    //let calcNrOfBins (lowerScanLimit,upperScanLimit)  = 
    //    ((lowerScanLimit-upperScanLimit) / 100.)
    //    |> abs
    //    |> ceil
    //    |> int

    /////
    //let spectrumToBinnedArray (qMostAbundandPepsMin, qMostAbundandPepsMax) (lowerScanLimit,upperScanLimit) (spectrumData:(float*float) []) = //(scanlimits:float*float) 
    //    /// TODO cut off specData < lowerlimit und specdata > upperlimit
    //    if spectrumData.Length = 0 then [||]
    //    else
    //    let spectrumData = spectrumData |> Array.filter (fun (mz,intensity) -> mz >= lowerScanLimit && mz <= upperScanLimit) 
    //    // Number of bins with width of 100 Thomson (th = m/z)
    //    let nrOfBins = 
    //        calcNrOfBins (lowerScanLimit,upperScanLimit)
    //    ///
    //    let binnedData = 
    //        let tmp = Array.create nrOfBins []
    //        for i = 0 to spectrumData.Length-1 do
    //            let binNumber = 
    //                fst spectrumData.[i] / 100.
    //                |> int 
    //            tmp.[binNumber-1] <- spectrumData.[i] :: tmp.[binNumber-1]
    //        tmp
    //    ///
    //    let intensitySortedBins = 
    //        binnedData
    //        |> Array.map (fun binL -> 
    //                        binL
    //                        |> List.sortByDescending (fun (mzData,intensity) -> intensity)
    //                        |> (fun x -> if x.Length >= qMostAbundandPepsMax then
    //                                          List.take qMostAbundandPepsMax x
    //                                     else List.take x.Length x
    //                           )
    //                     ) 
    //    ///
    //    let jaggedIntensitySortedBins = 
    //        [|for q = qMostAbundandPepsMin to qMostAbundandPepsMax do
    //            yield
    //               q, [for l = 0 to nrOfBins-1 do
    //                    let currentbin = intensitySortedBins.[l]
                    
    //                    if currentbin.Length > 0 && currentbin.Length >= q then              
    //                        yield currentbin 
    //                              |> List.take q
    //                              |> List.sortBy (fun (mzData,intensity) -> mzData)
                             
    //                    else
    //                        yield currentbin 
    //                              |> List.sortBy (fun (mzData,intensity) -> mzData)     
    //                ]
        
    //        |]
    //    jaggedIntensitySortedBins
      
    /////
    //let predictOf (lowerScanLimit,upperScanLimit) chargeState (fragments:PeakFamily<TaggedMass.TaggedMass> list)  =
    //    let predictPeak charge (taggedMass: TaggedMass.TaggedMass) = 
    //        TaggedPeak.TaggedPeak(taggedMass.Iontype, (Mass.toMZ taggedMass.Mass charge), nan)
    //    let computePeakFamily charge fragments = 
    //        let mainPeak = predictPeak charge fragments.MainPeak
    //        let dependentPeaks =
    //            fragments.DependentPeaks
    //            |> List.fold (fun acc dependent -> 
    //                             if charge <= 1. then 
    //                                 predictPeak charge dependent 
    //                                 :: acc
    //                             else 
    //                                 acc
    //                          ) []
    //        createPeakFamily mainPeak dependentPeaks
    //    let rec recloop ions charge (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
    //        match fragments with
    //        | fragments::rest ->  
    //            if charge <= 1. then 
    //                let tempIons = computePeakFamily 1. fragments 
    //                recloop (tempIons::ions) charge rest        
    //            else 

    //            let tempIons = 
    //                [

    //                for z = 1 to 2 do 
    //                    yield computePeakFamily (float z) fragments
    //                ]
    //            recloop (tempIons@ions) charge rest
    //        | [] -> ions       
    //    recloop [] chargeState fragments
    //    |> List.sortBy (fun peak -> peak.MainPeak.Mz)
    //    |> List.toArray
                        

    /////
    //let findBinIdx mz =
    //    (mz / 100. |> int) - 1

    //let hasMatchingPeakMZTol (mzMatchingTolerance: float) (tarmz: float) (observBinnedSpect: ((float*float) list) list)  =
    //    ///
    //    let targetBin = findBinIdx tarmz
    //    ///
    //    match List.tryFind ( fun (mz,int) -> ( (mz-tarmz) |> abs ) <= mzMatchingTolerance) observBinnedSpect.[targetBin]  with
    //    | Some _    -> true
    //    | None      -> false

    //let hasMatchingPeakPPM (mzMatchingTolerancePPM: float) (tarmz: float) (observBinnedSpect: ((float*float) list) list)  =
    //    ///
    //    let targetBin = findBinIdx tarmz
    //    if targetBin < 0 then false 
    //    else
    //    ///
    //    match List.tryFind ( fun (mz,int) -> ( (mz-tarmz) |> abs ) <= Mass.deltaMassByPpm mzMatchingTolerancePPM tarmz) observBinnedSpect.[targetBin]  with
    //    | Some _    -> true
    //    | None      -> false

    /////
    //let countMatches (scanlimits:float*float) (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (observBinnedSpect:int* ((float*float) list) list ) =    
    //    ///
    //    let lowerLimit,upperLimit = scanlimits
    //    let (q, binnedSpec) = observBinnedSpect

    //    /// declare variables
    //    let mutable nWithOutNeutralAcc  = 0
    //    let mutable nWithNeutralAcc     = 0
    //    let mutable kWithOutNeutralAcc  = 0
    //    let mutable kWithNeutralAcc     = 0
    //    let matchAcc = ResizeArray()
    //    let rec findMatches (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (observBinnedSpect: ((float*float) list) list ) idx nWithOutNeutral nWithNeutral kWithOutNeutral kWithNeutral matchesWithOutNeutral =
    //        if idx = theoSpect.Length then
    //            nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral, matchesWithOutNeutral
    //        else
    //            let currentPeakFam = theoSpect.[idx]
    //            //
    //            nWithOutNeutralAcc  <- 0
    //            nWithNeutralAcc     <- 0
    //            //
    //            kWithOutNeutralAcc  <- 0
    //            kWithNeutralAcc     <- 0

    //            if currentPeakFam.MainPeak.Mz >= lowerLimit && currentPeakFam.MainPeak.Mz <= upperLimit then                
    //                if hasMatchingPeakPPM mzMatchingTolerance currentPeakFam.MainPeak.Mz observBinnedSpect = true then
    //                    // raise ns
    //                    nWithOutNeutralAcc  <- 1
    //                    nWithNeutralAcc     <- 1     
    //                    // raise ks                          
    //                    kWithOutNeutralAcc  <- 1
    //                    kWithNeutralAcc     <- 1
    //                    matchAcc.Add currentPeakFam.MainPeak
    //                    currentPeakFam.DependentPeaks 
    //                    |> List.iter ( fun dependentAnnotPeak ->
    //                                    // Peak not within bounds loop further. 
    //                                    if dependentAnnotPeak.Mz >= lowerLimit && dependentAnnotPeak.Mz <= upperLimit then 
    //                                         if hasMatchingPeakPPM mzMatchingTolerance dependentAnnotPeak.Mz observBinnedSpect = true then
    //                                         //Successful match: raise n and k
    //                                            match dependentAnnotPeak.Iontype with
    //                                            //If NeutralLoss occurs only raise the neutralLossAccs 
    //                                            | iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
    //                                                nWithNeutralAcc     <- nWithNeutralAcc    + 1                               
    //                                                kWithNeutralAcc     <- kWithNeutralAcc    + 1
    //                                            | _                      ->
    //                                                // raise ns
    //                                                nWithOutNeutralAcc  <- nWithOutNeutralAcc + 1
    //                                                nWithNeutralAcc     <- nWithNeutralAcc    + 1     
    //                                                // raise ks                      
    //                                                kWithOutNeutralAcc  <- kWithOutNeutralAcc + 1
    //                                                kWithNeutralAcc     <- kWithNeutralAcc    + 1       
    //                                                matchAcc.Add dependentAnnotPeak
    //                                         else 
    //                                         //Unsuccessful match: raise n
    //                                            match dependentAnnotPeak.Iontype with 
    //                                            |  iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
    //                                                nWithNeutralAcc     <- nWithNeutralAcc    + 1

    //                                            | _                       ->      
    //                                                nWithOutNeutralAcc  <- nWithOutNeutralAcc + 1
    //                                                nWithNeutralAcc     <- nWithNeutralAcc    + 1   
    //                                    else 
    //                                        nWithOutNeutralAcc  <- nWithOutNeutralAcc + 1
    //                                        nWithNeutralAcc     <- nWithNeutralAcc    + 1  
    //                                  )  
    //                    findMatches mzMatchingTolerance theoSpect observBinnedSpect (idx+1) (nWithOutNeutral + nWithOutNeutralAcc) (nWithNeutral + nWithNeutralAcc) (kWithOutNeutral + kWithOutNeutralAcc) (kWithNeutralAcc + kWithNeutral) matchesWithOutNeutral
    //                else 
    //                    // No matching Peak. Raise n.
    //                    findMatches mzMatchingTolerance theoSpect observBinnedSpect (idx+1) (nWithOutNeutral+1) (nWithNeutral+1) kWithOutNeutral kWithNeutral matchAcc
    //            else 
    //                // Peak not within bounds raise n 
    //                findMatches mzMatchingTolerance theoSpect observBinnedSpect (idx+1) (nWithOutNeutral+1) (nWithNeutral+1) kWithOutNeutral kWithNeutral matchAcc

    //    let (nWithOutNeutral, nWithNeutral, kWithOutNeutral,kWithNeutral,matchesWithOutNeutral)  = findMatches mzMatchingTolerance theoSpect binnedSpec 0 0 0 0 0 matchAcc
    //    (nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral,matchesWithOutNeutral) 

    //let scoreFuncImpl n k topx =
    //    let topx' = float topx
    //    let mutable acc = 0.
    //    for j = k to n do 
    //        let binomialCoeff = MathNet.Numerics.SpecialFunctions.Binomial(n,j)
    //        let p1 = (topx'/100.)**(j |> float)
    //        let p2 = (1.-(topx'/100.)) ** (float (n-j))
    //        acc <- acc + (binomialCoeff * p1 * p2)
    //    -10.*Math.Log10(acc) 
    
    //let massCorrection mass = 
    //    0.024*(mass - 600.)

    //let modCorrection nMod = 
    //    match nMod with 
    //    | 0  -> 42.
    //    | 1  -> 28.
    //    | 2  -> 22.
    //    | 3  -> 16.
    //    | 4  -> 9.
    //    | 5  -> 5.
    //    | 6  -> 2.
    //    | _  -> 0.
    
    //let cleavageCorrection nCleave consecutiveCleavage =
    //    match nCleave with 
    //    | 0 -> 53.2
    //    | 1 -> 
    //        if consecutiveCleavage then 
    //            42.1
    //        else 
    //            31.1
    //    | 2 -> 17. 
    //    |_  -> 0.

        
    //let scoreTheoVsRecordedSpec scanlimits matchingTolPPM (lookUpResult:LookUpResult<'a>) theoSpec binnedRecSpec= 
    //    binnedRecSpec
    //    |> Array.map (fun (q,l) -> q, countMatches scanlimits matchingTolPPM theoSpec (q,l))
    //    |> Array.map (fun (q,(nWo,nWi,kWo,kWi,matchesWithoutNeutral)) -> 
                                
    //                            let massCorr = massCorrection lookUpResult.Mass
    //                            let modCorr  = modCorrection 0
    //                            let cleavageCorr = cleavageCorrection 0 false
    //                            let rawScore1 = scoreFuncImpl (nWo) kWo q
    //                            let rawScore2 = scoreFuncImpl (nWi) kWi q 
    //                            if rawScore1 > rawScore2 then 
    //                                createMatchingScore (rawScore1 + massCorr + modCorr + cleavageCorr - 100.) nWo kWo q (matchesWithoutNeutral |> Seq.toList)
    //                            else 
    //                                createMatchingScore (rawScore2 + massCorr + modCorr + cleavageCorr - 100.) nWi kWi q (matchesWithoutNeutral |> Seq.toList)
    //                    )
    //    |> Array.max 
    
    /////
    //let calcDeltaAndorScore (sourceList:SearchEngineResult<MatchingScore> list) =
    //    match sourceList with
    //    | h1::rest -> 
    //        sourceList
    //        |> List.map
    //            (fun sls ->
    //                let deltaAndroScore = (h1.Score.Score - sls.Score.Score) / h1.Score.Score
    //                { sls with DeltaCN = deltaAndroScore } )      
    //    | []       -> []
    
    /////
    //let calcDeltaAndorScoreBy (sourceList:SearchEngineResult<float> list) =
    //    match sourceList with
    //    | h1::rest -> 
    //        sourceList
    //        |> List.map
    //            (fun sls ->
    //                let deltaAndroScore = (h1.Score - sls.Score) / h1.Score
    //                { sls with DeltaCN = deltaAndroScore } )      
    //    | []       -> []
    
    /////
    //let calcAndromedaLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =
    //    let fCharge = float chargeState
    //    //
    //    let binnedRecSpec = 
    //        spectrum 
    //        |> Array.map (fun peak -> peak.Mz,peak.Intensity)
    //        |> spectrumToBinnedArray qMinAndMax (scanlimits:float*float) 
    //    //
    //    let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
    //    //
    //    let ides = 
    //        possiblePeptideInfos
    //        |> List.fold (fun acc lookUpResult -> 
    //                            let sequence = lookUpResult.BioSequence
    //                            let seqL = sequence.Length
    //                            let ionSeries = calcIonSeries massfunction sequence
    //                            // 
    //                            let theoSpec = 
    //                                predictOf scanlimits fCharge ionSeries
    //                            // Calculates score for real sequence
    //                            let targetScore = scoreTheoVsRecordedSpec scanlimits matchingTolPPM lookUpResult theoSpec binnedRecSpec
    //                            let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore nan

    //                            // Calculates score for reversed decoy sequence
    //                            let sequence_decoy = sequence |> List.rev
    //                            let ionSeries_decoy = calcIonSeries massfunction sequence_decoy
    //                            let theoSpecDecoy = predictOf scanlimits fCharge ionSeries_decoy
    //                            let decoyScore =  scoreTheoVsRecordedSpec scanlimits matchingTolPPM lookUpResult theoSpecDecoy binnedRecSpec
    //                            let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore nan 
    //                            targetResult :: decoyResult :: acc                        
    //                     ) []
    //    calcDeltaAndorScore (ides  |> List.sortBy (fun sls -> - sls.Score.Score))  

    /////
    //let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
    //    TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos
  
    /////
    //let inline calcAndromedaScore qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
    //    let fCharge = float chargeState
    //    //
    //    let binnedRecSpec = 
    //        spectrum 
    //        |> Array.map (fun peak -> peak.Mz,peak.Intensity)
    //        |> spectrumToBinnedArray qMinAndMax scanlimits
    //    //
    //    let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
    //    let ides = 
    //        theoreticalSpectra 
    //        |> List.fold (fun acc theoreticalSpectrum -> 
    //                            let lookUpResult = theoreticalSpectrum.LookUpResult
    //                            let seqL = lookUpResult.BioSequence.Length
 
    //                            // Calculates score for target sequence
    //                            let theoSpec = theoreticalSpectrum.TheoSpec
    //                            let targetScore = scoreTheoVsRecordedSpec scanlimits matchingTolPPM lookUpResult theoSpec binnedRecSpec
    //                            let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan
                                
    //                            let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
    //                            let decoyScore =  scoreTheoVsRecordedSpec scanlimits matchingTolPPM lookUpResult theoSpec_Decoy binnedRecSpec
    //                            let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan
    //                            targetResult :: decoyResult :: acc     
    //                 ) []
    //    calcDeltaAndorScoreBy (ides  |> List.sortBy (fun sls -> - sls.Score))       
    