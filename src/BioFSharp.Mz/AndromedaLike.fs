namespace BioFSharp.Mz

open System
open BioFSharp
open BioFSharp.Mz.Peaks
open SearchDB
open Fragmentation
open TheoreticalSpectra
open SearchEngineResult
open FSharpAux.Array
//open MathNet.Numerics

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

    /// Counts how many peaks present in spectrum data within a a given mzWindow (mzwindow/2 +/- ) have a higher 
    /// intensity then the peak at position spectrumData.[peakIdx]. 
    let private moreIntensePeaksWithin mzWindow (spectrumData:(float*float) []) peakIdx =
        let halfWindowSize = mzWindow / 2. 
        let currentPeakMz,currentPeakIntensity = spectrumData.[peakIdx] 
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
            loop peakIdx 0 
        let leftSidePeaks = 
            let rec loop currentIdx moreIntensePeaks = 
                if currentIdx < 0 then moreIntensePeaks
                else 
                    let tempPeakMz,tempPeakInt = spectrumData.[currentIdx] 
                    if tempPeakMz < lowerBorder then moreIntensePeaks
                    elif tempPeakInt > currentPeakIntensity then 
                         loop (currentIdx-1) (moreIntensePeaks+1) 
                    else loop (currentIdx-1) moreIntensePeaks 
            loop peakIdx 0         
        createRatedPeak (currentPeakMz,currentPeakIntensity) (leftSidePeaks+rightSidePeaks) 
    
    /// Applies the function "moreIntensePeaksWithin" to every peak in spectrumData, thus rating every peak by counting how many more abundant peaks
    /// are present within the given mzWindow (mzwindow/2 +/- ). Filters out peaks that are not within the scanLimits and show to have more abundand peaks than 
    /// defined by qMostAbundandpepsMax in their neighbourhood.
    let private ratedSpectrum mzWindow (qMostAbundandPepsMin, qMostAbundandPepsMax) (lowerScanLimit,upperScanLimit) (spectrumData:(float*float) []) = 
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
          
    /// Converts the fragmentIon ladder to a theoretical spectrum at a given charge state. Filters out all theoretical peaks, that lie
    /// outside the given lower and upperScanLimits.
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

    /// Searches for a Peak in ratedSpectrum that matches the target mz (tarmz) within a certain tolerance defined by mzmatchingTolerance.
    let private hasMatchingPeakMZTol (mzMatchingTolerance: float) (tarmz: float) (ratedSpectrum: RatedPeak []) =        
        match binarySearchIndexBy  (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= mzMatchingTolerance then 0 
                                elif tarmz < mz then -1
                                else 1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | x             ->  None 

    /// Searches for a Peak in ratedSpectrum that matches the target mz (tarmz) within a certain tolerance expressed in parts per million of the tarmz, defined by mzMatchingTolerancePPM.
    let private hasMatchingPeakPPM (mzMatchingTolerancePPM: float) (tarmz: float) (ratedSpectrum: RatedPeak []) = 
        match binarySearchIndexBy (fun (rPeak: RatedPeak) ->
                                let mz = fst rPeak.Peak  
                                if ((mz-tarmz) |> abs) <= Mass.deltaMassByPpm mzMatchingTolerancePPM tarmz then 0 
                                elif tarmz < mz then -1
                                else 1
                            ) ratedSpectrum with
        | x when x >= 0 -> Some ratedSpectrum.[x] 
        | _             -> None 

    /// Iterates theoSpect and counts the number of common peaks in ratedspectrum. A peak is considered as common if the mz value lies within a user defined matching tolerance. This procedure
    /// is carried out (qMostAbundandPepsMax-qMostAbundandPepsMin+1) times every time with a different qThreshold to account for peak intensity. A peak is counted as matching if it has
    /// less more abundant peaks in the neighbourhood than the current qThreshold.
    let private countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) (mzMatchingTolerance: float) (theoSpect: PeakFamily<TaggedPeak.TaggedPeak> []) (ratedSpectrum: RatedPeak []) =    
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


    /// Helper function to calculate the AndromedaScore
    let private lnProb n k lnp lnq = 
        let s1 = -k * lnp - (n - k) * lnq - FSharp.Stats.SpecialFunctions.Gamma.gammaLn(n + 1.) 
                    + FSharp.Stats.SpecialFunctions.Gamma.gammaLn(k + 1.) + FSharp.Stats.SpecialFunctions.Gamma.gammaLn(n - k + 1.) 
        s1 
    
    /// Calculates the andromedaLike Score based on offered peaks n and matched peaks k, as well as the q value threshold topx
    let private scoreFuncImpl n k topx = 
        let fTopx = float topx 
        let p1 = Math.Min(fTopx/100.0,0.5)
        let lnp = Math.Log(p1)
        let lnq = Math.Log(1.0-p1)
        let mutable acc = 0.
        for i = k to n do 
            acc <- acc + Math.Exp(-lnProb (float n) (float i) lnp lnq)
        acc <- -Math.Log(acc)
        10. * acc / Math.Log(10.)

    /// Correction factor applided to the andromeda score based on the observed precursor Mz.
    /// Values are taken from the original andromeda release.
    let private mzCorrection mass = 
        0.024*(mass - 600.)

    /// Correction factor applided to the andromeda score based on the amount modifications.
    /// Values are taken from the original andromeda release.
    let private modCorrection nMod = 
        match nMod with 
        | 0  -> 42.
        | 1  -> 28.
        | 2  -> 22.
        | 3  -> 16.
        | 4  -> 9.
        | 5  -> 5.
        | 6  -> 2.
        | _  -> 0.
    
    /// Correction factor applided to the andromeda score based on the amount of missCleavages.
    /// Values are taken from the original andromeda release.
    let private cleavageCorrection nCleave consecutiveCleavage =
        match nCleave with 
        | 0 -> 53.2
        | 1 -> 
            if consecutiveCleavage then 
                42.1
            else 
                31.1
        | 2 -> 17. 
        |_  -> 0.


    /// Computes the andromedalike score of a theoretical spectrum vs. the ratedSpectrum. Subsequently applies mzCorrection - modCorrection and cleavageCorrection to the score. 
    let private scoreTheoVsRecordedSpec (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax)  matchingTolPPM  precursorMZ (*charge*) theoSpec (ratedSpectrum: RatedPeak []) =     
        countMatches (lowerScanLimit,upperScanLimit) (qMostAbundandPepsMin, qMostAbundandPepsMax) matchingTolPPM theoSpec ratedSpectrum
        |> Array.map (fun (q,(nWo,nWi,kWo,kWi)) -> 
                        let massCorr = mzCorrection precursorMZ //(Mass.toMZ lookUpResult.Mass charge)
                        let modCorr  = modCorrection 0
                        let cleavageCorr = cleavageCorrection 0 false
                        let rawScore1 = scoreFuncImpl (nWo) kWo q
                        let rawScore2 = scoreFuncImpl (nWi) kWi q 

                        if rawScore1 > rawScore2 then 
                            let finalScore = 
                                let tmp = (rawScore1 + massCorr + modCorr + cleavageCorr - 100.)
                                if tmp < 0. then 0. else tmp 
                            createMatchingScore finalScore nWo kWo q
                        else 
                            let finalScore = 
                                let tmp = (rawScore2 + massCorr + modCorr + cleavageCorr - 100.)
                                if tmp < 0. then 0. else tmp 
                            createMatchingScore finalScore nWi kWi q 
                    )
        |> Array.maxBy (fun x -> x.Score) 
                            

    /// Converts the fragment ion ladders to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos
 
    /// Calculates the AndromedaLike Scores for all peptides (possiblePeptideInfos).
    let calcAndromedaLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =
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
                                let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec ratedSpec
                                let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan nan

                                // Calculates score for reversed decoy sequence
                                let sequence_decoy = sequence |> List.rev
                                let ionSeries_decoy = calcIonSeries massfunction sequence_decoy
                                let theoSpecDecoy = predictOf scanlimits fCharge ionSeries_decoy
                                let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpecDecoy ratedSpec
                                let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan nan
                                targetResult :: decoyResult :: acc                        
                         ) []
        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))  
        |> calcNormDeltaNext


    /// Calculates the AndromedaLike scores for all theoretical spectra.
    let calcAndromedaScore qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
        
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
                                let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec binnedRecSpec
                                let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan nan
                                
                                let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
                                let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec_Decoy binnedRecSpec
                                let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan nan
                                targetResult :: decoyResult :: acc     
                     ) []
        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))       
        |> calcNormDeltaNext

   
    /// Calculates the AndromedaLike scores for all theoretical spectra. Implemented using Async parallel. 
    let calcAndromedaScoreParallel qMinAndMax scanlimits matchingTolPPM (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<PeakFamily<TaggedPeak.TaggedPeak>[]>> ) spectrumID =
        
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
            let targetScore = scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec binnedRecSpec
            let targetResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL targetScore.Score nan nan
            
            let theoSpec_Decoy = theoreticalSpectrum.DecoyTheoSpec
            let decoyScore =  scoreTheoVsRecordedSpec scanlimits qMinAndMax matchingTolPPM isolationWindowTargetMz theoSpec_Decoy binnedRecSpec
            let decoyResult = createSearchEngineResult SearchEngineResult.SearchEngine.AndromedaLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass seqL decoyScore.Score nan nan
            [targetResult;decoyResult]
        
        let ides = 
            [for i in theoreticalSpectra ->  async {return computeScores i}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> List.concat    
        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))       
        |> calcNormDeltaNext

