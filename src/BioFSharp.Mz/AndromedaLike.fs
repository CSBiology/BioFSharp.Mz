namespace BioFSharp.Mz

open System
open BioFSharp
open BioFSharp.Mz.Peaks
open SearchDB
open Fragmentation
open TheoreticalSpectra
open SearchEngineResult

open MathNet.Numerics

module AndromedaLike =
    
    type MatchingScore = {
        Score   : float 
        N       : int 
        K       : int
        Q       : int
        //Matches : PeakAnnotation list 
        }
    
    let createMatchingScore score n k q = {
     Score = score; N = n ; K = k; Q = q ;//Matches = matches
    }

    [<Struct>]
    type RatedPeak = {
        Peak: float*float
        MoreIntensePeaksWithinWindow: int
        }

    let createRatedPeak peak moreIntensePeaks = {
        Peak = peak; MoreIntensePeaksWithinWindow= moreIntensePeaks }


    let moreIntensePeaksWithin mzWindow (spectrumData:(float*float) []) currentIdx =
        let halfWindowSize = mzWindow / 2. 
        let currentPeakMz,currentPeakIntensity = spectrumData.[currentIdx] 
        let upperBorder = currentPeakMz + halfWindowSize
        let lowerBorder = currentPeakMz - halfWindowSize
        let rightSidePeaks =
            let rec loop currentIdx moreIntensePeaks = 
                if currentIdx = spectrumData.Length then moreIntensePeaks
                else 
                    let tempPeakMz,tempPeakInt = spectrumData.[currentIdx] 
                    if tempPeakMz > upperBorder then moreIntensePeaks
                    elif tempPeakInt > currentPeakIntensity then 
                         loop (currentIdx+1) (moreIntensePeaks+1) 
                    else loop (currentIdx+1) moreIntensePeaks 
            loop currentIdx 0 
        let leftSidePeaks = 
            let rec loop currentIdx moreIntensePeaks = 
                if currentIdx < 0 then moreIntensePeaks
                else 
                    let tempPeakMz,tempPeakInt = spectrumData.[currentIdx] 
                    if tempPeakMz < lowerBorder then moreIntensePeaks
                    elif tempPeakInt > currentPeakIntensity then 
                         loop (currentIdx-1) (moreIntensePeaks+1) 
                    else loop (currentIdx-1) moreIntensePeaks 
            loop currentIdx 0         
        createRatedPeak (currentPeakMz,currentPeakIntensity) (leftSidePeaks+rightSidePeaks) 
    
    ///
    let ratedSpectrum mzWindow (qMostAbundandPepsMin, qMostAbundandPepsMax) (lowerScanLimit,upperScanLimit) (spectrumData:(float*float) []) = 
        if spectrumData.Length = 0 then [||]
        else
        let filteredAndRanked = 
            spectrumData
            |> Array.mapi   (fun i (mzData,intensityData) -> 
                                moreIntensePeaksWithin mzWindow spectrumData i      
                            )   
            |> Array.filter (fun rPeak -> 
                                let mzData = fst rPeak.Peak  
                                // scanLimits vllt raus
                                mzData >= lowerScanLimit && mzData <= upperScanLimit && rPeak.MoreIntensePeaksWithinWindow < qMostAbundandPepsMax
                            )                        
        filteredAndRanked
          
    ///
    let predictOf (lowerScanLimit,upperScanLimit) chargeState (fragments:PeakFamily<TaggedMass.TaggedMass> list)  =
        let predictPeak charge (taggedMass: TaggedMass.TaggedMass) = 
            TaggedPeak.TaggedPeak(taggedMass.Iontype, (Mass.toMZ taggedMass.Mass charge), nan)
        let computePeakFamily charge fragments = 
            let mainPeak = predictPeak charge fragments.MainPeak
            let dependentPeaks =
                fragments.DependentPeaks
                |> List.fold (fun acc dependent -> 
                                 if charge <= 1. then 
                                     predictPeak charge dependent 
                                     :: acc
                                 else 
                                     acc
                              ) []
            createPeakFamily mainPeak dependentPeaks
        let rec recloop ions charge (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
            match fragments with
            | fragments::rest ->  
                if charge <= 1. then 
                    let tempIons = computePeakFamily 1. fragments 
                    recloop (tempIons::ions) charge rest        
                else 

                let tempIons = 
                    [
                    for z = 1 to 2 do 
                        yield computePeakFamily (float z) fragments
                    ]
                recloop (tempIons@ions) charge rest
            | [] -> ions       
        recloop [] chargeState fragments
        //|> List.sortBy (fun peak -> peak.MainPeak.Mz)
        |> List.toArray
                        
    ///
    let binarySearch f (arr: 'a []) = 
        if arr.Length = 0 then
            0
        else
            let rec loop lower upper = 
                if lower > upper then ~~~ lower 
                else
                    let middle = lower + ((upper - lower) / 2)
                    let comparisonResult = f arr.[middle]   
                    if comparisonResult = 0 then
                        middle
                    elif comparisonResult < 0 then
                        loop lower (middle - 1)
                    else
                        loop (middle + 1) upper           
            loop 0 (arr.Length-1) 

    ///
    let hasMatchingPeakMZTol (mzMatchingTolerance: float) (tarmz: float) (ratedSpectrum: RatedPeak []) =        
        match binarySearch  (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= mzMatchingTolerance then 0 
                                elif mz < tarmz then 1
                                else -1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | x            ->  None 

    ///
    let hasMatchingPeakPPM (mzMatchingTolerancePPM: float) (tarmz: float) (ratedSpectrum: RatedPeak []) = 
        match binarySearch  (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= Mass.deltaMassByPpm mzMatchingTolerancePPM tarmz then 0 
                                elif mz < tarmz then 1
                                else -1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | _            -> None 

    ///
    let countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (ratedSpectrum: RatedPeak []) =    
        let countArray = Array.init (qMostAbundandPepsMax-qMostAbundandPepsMin+1) (fun i -> qMostAbundandPepsMin+i,(0,0,0,0))
        let findMaxIterator (rPeak:RatedPeak) = 
            match (1+rPeak.MoreIntensePeaksWithinWindow)-qMostAbundandPepsMin with 
            | x when x <= 0 -> 0
            | x when x >= countArray.Length-1 -> x
            | x -> x         
        let rec findMatches (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (ratedSpectrum: RatedPeak []) idx =
            if idx = theoSpect.Length then
                countArray
            else
                let currentPeakFam = theoSpect.[idx]
                if currentPeakFam.MainPeak.Mz >= lowerScanLimit && currentPeakFam.MainPeak.Mz <= upperScanLimit then                
                    match hasMatchingPeakPPM mzMatchingTolerance currentPeakFam.MainPeak.Mz ratedSpectrum with 
                    | Some rPeak ->
                        let maxIterator = findMaxIterator rPeak
                        //raise ns and ks 
                        for i = 0 to countArray.Length-1 do
                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                            if i >= maxIterator then
                                countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral+1, kWithNeutral+1))      
                            else 
                                countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral, kWithNeutral))
                        currentPeakFam.DependentPeaks 
                        |> List.iter ( fun dependentAnnotPeak ->
                                        // Peak not within bounds loop further. 
                                        if dependentAnnotPeak.Mz >= lowerScanLimit && dependentAnnotPeak.Mz <= upperScanLimit then 
                                             match hasMatchingPeakPPM mzMatchingTolerance dependentAnnotPeak.Mz ratedSpectrum with 
                                             //Successful match: raise n and k
                                             | Some rPeak -> 
                                                 match dependentAnnotPeak.Iontype with
                                                    //If NeutralLoss occurs only raise the neutralLossAccs 
                                                    | iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
                                                        let maxIterator = findMaxIterator rPeak
                                                        //raise ns and ks 
                                                        for i = 0 to countArray.Length-1 do
                                                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                                                            if i >= maxIterator then
                                                                countArray.[i] <- (q,(nWithOutNeutral, nWithNeutral+1, kWithOutNeutral, kWithNeutral+1))      
                                                            else 
                                                                countArray.[i] <- (q,(nWithOutNeutral, nWithNeutral+1, kWithOutNeutral, kWithNeutral))
                                                    | _                      ->
                                                        //raise ns and ks 
                                                        for i = 0 to countArray.Length-1 do
                                                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                                                            if i >= maxIterator then
                                                                countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral+1, kWithNeutral+1))      
                                                            else 
                                                                countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral, kWithNeutral))       
                                             | None      -> 
                                             //Unsuccessful match: raise n
                                                match dependentAnnotPeak.Iontype with 
                                                |  iontype when iontype.HasFlag(Ions.IonTypeFlag.Neutral) ->
                                                        for i = 0 to countArray.Length-1 do
                                                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                                                            countArray.[i] <- (q,(nWithOutNeutral, nWithNeutral+1, kWithOutNeutral, kWithNeutral))  
                                                | _                       ->      
                                                        for i = 0 to countArray.Length-1 do
                                                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                                                            countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral, kWithNeutral))   
                                        )  
                        findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1) 
                    | None ->
                        // No matching Peak. Raise ns. 
                        for i = 0 to countArray.Length-1 do
                            let (q,(nWithOutNeutral, nWithNeutral, kWithOutNeutral, kWithNeutral)) = countArray.[i] 
                            countArray.[i] <- (q,(nWithOutNeutral+1, nWithNeutral+1, kWithOutNeutral, kWithNeutral))
                        findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1)
                else 
                    // Peak not within bounds raise n 
                    findMatches mzMatchingTolerance theoSpect ratedSpectrum (idx+1)

        findMatches mzMatchingTolerance theoSpect ratedSpectrum 0


    ///
    let lnProb n k lnp lnq = 
        let s1 = -k * lnp - (n - k) * lnq - MathNet.Numerics.SpecialFunctions.GammaLn(n + 1.) 
                    + MathNet.Numerics.SpecialFunctions.GammaLn(k + 1.) + MathNet.Numerics.SpecialFunctions.GammaLn(n - k + 1.) 
        s1 
    
    ///
    let log10' = Math.Log(10.)

    ///
    let scoreFuncImpl n k topx = 
        let fTopx = float topx 
        let p1 = Math.Min(fTopx/100.0,0.5)
        let lnp = Math.Log(p1)
        let lnq = Math.Log(1.0-p1)
        let mutable acc = 0.
        for i = k to n do 
            acc <- acc + Math.Exp(-lnProb (float n) (float k) lnp lnq)
        acc <- -Math.Log(acc)
        10. * acc / log10'

    //let scoreFuncImpl n k topx =
    //    let topx' = float topx
    //    let mutable acc = 0.
    //    for j = k to n do 
    //        let binomialCoeff = MathNet.Numerics.SpecialFunctions.Binomial(n,j)
    //        let p1 = (topx'/100.)**(j |> float)
    //        let p2 = (1.-(topx'/100.)) ** (float (n-j))
    //        acc <- acc + (binomialCoeff * p1 * p2)
    //    -10.*Math.Log10(acc) 
    
    let massCorrection mass = 
        0.024*(mass - 600.)

    let modCorrection nMod = 
        match nMod with 
        | 0  -> 42.
        | 1  -> 28.
        | 2  -> 22.
        | 3  -> 16.
        | 4  -> 9.
        | 5  -> 5.
        | 6  -> 2.
        | _  -> 0.
    
    let cleavageCorrection nCleave consecutiveCleavage =
        match nCleave with 
        | 0 -> 53.2
        | 1 -> 
            if consecutiveCleavage then 
                42.1
            else 
                31.1
        | 2 -> 17. 
        |_  -> 0.

        
    let scoreTheoVsRecordedSpec (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax)  matchingTolPPM  charge (lookUpResult:LookUpResult<'a>) theoSpec (ratedSpectrum: RatedPeak []) =     
        countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) matchingTolPPM theoSpec ratedSpectrum
        |> Array.map (fun (q,(nWo,nWi,kWo,kWi)) -> 
                        let massCorr = massCorrection lookUpResult.Mass //(Mass.toMZ lookUpResult.Mass charge)
                        let modCorr  = modCorrection 0
                        let cleavageCorr = cleavageCorrection 0 false
                        let rawScore1 = scoreFuncImpl (nWo) kWo q
                        let rawScore2 = scoreFuncImpl (nWi) kWi q 
                        if rawScore1 > rawScore2 then 
                            createMatchingScore (rawScore1 + massCorr + modCorr + cleavageCorr - 100.) nWo kWo q
                        else 
                            createMatchingScore (rawScore2 + massCorr + modCorr + cleavageCorr - 100.) nWi kWi q 
                    )
        |> Array.max 

    ///
    let calcDeltaAndorScore (sourceList:SearchEngineResult<MatchingScore> list) =
        match sourceList with
        | h1::rest -> 
            sourceList
            |> List.map
                (fun sls ->
                    let deltaAndroScore = (h1.Score.Score - sls.Score.Score) / h1.Score.Score
                    { sls with DeltaCN = deltaAndroScore } )      
        | []       -> []
    
    ///
    let calcDeltaAndorScoreBy (sourceList:SearchEngineResult<float> list) =
        match sourceList with
        | h1::rest -> 
            sourceList
            |> List.map
                (fun sls ->
                    let deltaAndroScore = (h1.Score - sls.Score) / h1.Score
                    { sls with DeltaCN = deltaAndroScore } )      
        | []       -> []
    
    ///
    let calcAndromedaLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =
        let fCharge = float chargeState
        //
        let ratedSpec = 
            spectrum 
            |> Array.map (fun peak -> peak.Mz,peak.Intensity)
            |> ratedSpectrum 100. qMinAndMax scanlimits
        //
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        //
        let ides = 
            possiblePeptideInfos
            |> List.fold (fun acc lookUpResult -> 
                                let sequence = lookUpResult.BioSequence
                                let seqL = sequence.Length
                                let ionSeries = calcIonSeries massfunction sequence
                                // 
                                let theoSpec = 
                                    predictOf scanlimits fCharge ionSeries
                                // Calculates score for real sequence
                                let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpec ratedSpec
                                let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore nan

                                // Calculates score for reversed decoy sequence
                                let sequence_decoy = sequence |> List.rev
                                let ionSeries_decoy = calcIonSeries massfunction sequence_decoy
                                let theoSpecDecoy = predictOf scanlimits fCharge ionSeries_decoy
                                let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpecDecoy ratedSpec
                                let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore nan 
                                targetResult :: decoyResult :: acc                        
                         ) []
        calcDeltaAndorScore (ides  |> List.sortBy (fun sls -> - sls.Score.Score))  

    ///
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos
 
    ///
    let calcAndromedaScore qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
        
        let fCharge = float chargeState
        //
        let binnedRecSpec = 
            spectrum 
            |> Array.map (fun peak -> peak.Mz,peak.Intensity)
            |> ratedSpectrum 100. qMinAndMax scanlimits
        //
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        let ides = 
            theoreticalSpectra 
            |> List.fold (fun acc theoreticalSpectrum -> 
                                let lookUpResult = theoreticalSpectrum.LookUpResult
                                let seqL = lookUpResult.BioSequence.Length
 
                                // Calculates score for target sequence
                                let theoSpec = theoreticalSpectrum.TheoSpec
                                let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpec binnedRecSpec
                                let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan
                                
                                let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
                                let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpec_Decoy binnedRecSpec
                                let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan
                                targetResult :: decoyResult :: acc     
                     ) []
        calcDeltaAndorScoreBy (ides  |> List.sortBy (fun sls -> - sls.Score))       
   
   
    ///
    let calcAndromedaScoreParallel qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
        
        let fCharge = float chargeState
        //
        let binnedRecSpec = 
            spectrum 
            |> Array.map (fun peak -> peak.Mz,peak.Intensity)
            |> ratedSpectrum 100. qMinAndMax scanlimits
        //
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        let computeScores theoreticalSpectrum =
            let lookUpResult = theoreticalSpectrum.LookUpResult
            let seqL = lookUpResult.BioSequence.Length

            // Calculates score for target sequence
            let theoSpec = theoreticalSpectrum.TheoSpec
            let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpec binnedRecSpec
            let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan
            
            let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
            let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM fCharge lookUpResult theoSpec_Decoy binnedRecSpec
            let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan
            [targetResult;decoyResult]
        
        let ides = 
            [for i in theoreticalSpectra ->  async {return computeScores i}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> List.concat    
        calcDeltaAndorScoreBy (ides  |> List.sortBy (fun sls -> - sls.Score))       
    
