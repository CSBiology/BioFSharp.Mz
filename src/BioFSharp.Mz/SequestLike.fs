namespace BioFSharp.Mz

open BioFSharp
open SearchDB
open Fragmentation
open TheoreticalSpectra
open FSharpAux
open SearchEngineResult
open Ions
open FSharp.Stats
//open MathNet.Numerics
//open MathNet.Numerics.LinearAlgebra.Double

module SequestLike =

    /// normalize the intensities within a window to maximum of the window
    /// Attention shortens the array  (cuts)
    let private windowNormalizeIntensities (intensities:Vector<float>) (numberOfWindows:int) =
        // finds max within range
        let rec findMax (array:Vector<float>) (cMax:float) (lowerLimit:int) (counter:int)  =
            if counter < lowerLimit then
                cMax
            else
                let nMax = max array.[counter] cMax
                findMax array nMax lowerLimit (counter - 1)
        // sqrt intensities and normalize to maximum within lower-upper limit
        let normToMaxSqrtInPlace (array:Vector<float>) (lowerLimit:int) (upperLimit:int)  =
            let nMax = sqrt (findMax array 0.0 lowerLimit upperLimit)
            for i = lowerLimit to upperLimit do
                if nMax > 0. then                
                    array.[i] <- (sqrt array.[i]) / nMax 
                else
                    array.[i] <- 0.
        let windowSize =  (intensities.Length / numberOfWindows)    
        let tmpIntensities =            
            Vector.init (windowSize*numberOfWindows) (fun i -> intensities.[i])            
        for i = 1 to numberOfWindows do
            //printfn "window: %i lower: %i counter: %i " i (windowSize * (i - 1)) (windowSize * i - 1)
            normToMaxSqrtInPlace tmpIntensities (windowSize * (i - 1)) (windowSize * i - 1) 
        tmpIntensities
    
    /// Predicts the intensity of a theoretical peak based on the given charge and iontype.
    let private predictIntensitySimpleModel (ionType:Ions.IonTypeFlag) (charge:float) =
        match ionType with 
        | ionT when hasFlag Ions.IonTypeFlag.Unknown    ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Diagnostic ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Neutral    ionT -> 0.6 / charge
        | ionT when hasFlag Ions.IonTypeFlag.Immonium   ionT -> 0.6 / charge
        | ionT when hasFlag Ions.IonTypeFlag.lossNH3    ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.lossH2O    ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.A          ionT -> 0.2 / charge   
        | ionT when hasFlag Ions.IonTypeFlag.B          ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.C          ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.X          ionT -> 0.2 / charge
        | ionT when hasFlag Ions.IonTypeFlag.Y          ionT -> 1.  / charge
        | ionT when hasFlag Ions.IonTypeFlag.Z          ionT -> 0.2 / charge
        | _                                                  -> 0.2 / charge

    let private setVectorInplace (vector:Vector<float>) charge maxIndex lowerScanLimit (taggedMass:TaggedMass.TaggedMass) = 
        let index = int(System.Math.Round (Mass.toMZ taggedMass.Mass charge) ) - lowerScanLimit 
        if index < maxIndex-1 && index > -1 then
            vector.[index] <- max vector.[index] (predictIntensitySimpleModel (taggedMass.Iontype) charge) 

    /// Converts the fragment ion ladder to a theoretical Sequestlike spectrum at a given charge state. 
    let theoSpecOf (lowerScanLimit,upperScanLimit) (maxcharge:float) (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
        let lowerScanLimit = int lowerScanLimit
        let upperScanLimit = int upperScanLimit
        let maxIndex = upperScanLimit - lowerScanLimit + 1        
        let vector = Vector.create (maxIndex-1) 0.
        fragments 
        |> List.iter (fun p ->  
            let peaks = p.MainPeak::p.DependentPeaks
            [1. .. maxcharge]  
            |> List.iter (fun ch -> peaks |> List.iter (setVectorInplace vector ch maxIndex lowerScanLimit) ) 
            )
        vector

    /// Computes the autocorrelation of the vector +/- plusMinusMaxDelay.
    let private autoCorrelation (plusMinusMaxDelay:int) (vector:Vector<float>) =
        let shifted (vector:Vector<float>) (tau:int) =
            vector
            |> Vector.mapi
                (fun i x ->
                    let index = i - tau
                    if (index < 0) || (index > vector.Length - 1) then 
                        0.
                    else
                        vector.[index] )
        let rec accumVector (accum) (state:int) (max:int) =
            if state = max then
                accum
            else
                accumVector (accum + (shifted vector state)) (state - 1) (max)
        let emtyVector = Vector.zero vector.Length
        let plus  = accumVector emtyVector (plusMinusMaxDelay) (1)
        let minus = accumVector emtyVector (-1) (-plusMinusMaxDelay)
        (plus + minus)
        |> Vector.map (fun x ->  x / (float plusMinusMaxDelay * 2.))
       
    /////
    //let private createShiftedMatrix (plusMinusMaxDelay:int) (array:float[]) =
    //    let colNumber = array.Length
    //    Array2D.init (plusMinusMaxDelay * 2 + 1) colNumber
    //        (fun i ii ->
    //            let ni = (i - plusMinusMaxDelay)
    //            let index = ii - ni
    //            if (index < 0) || (index > colNumber - 1) then 
    //                0.
    //            else
    //                array.[index] )

    /// Converts the fragment ion ladder to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let predictOf (lowerScanLimit,upperScanLimit) charge ionSeries =         
        theoSpecOf (lowerScanLimit,upperScanLimit) charge ionSeries      
        

//    /// Measured spectrum to sequest-like normalized intensity array
//    /// ! Uses 10 as number of windows for window normalization (like in original sequest algorithm)    
//    let spectrumToNormalizedIntensityArray (scanlimits:int*int) (spectrum:PeakArray<_>) =
//        let lowerScanLimit,upperScanLimit = scanlimits    
//        let si  = PeakArray.peaksToNearestUnitDaltonBin spectrum lowerScanLimit upperScanLimit
//        let nsi = windowNormalizeIntensities si 10 // |> Seq.toArray
//        nsi

    /// Measured spectrum to sequest-like normalized intensity array
    /// minus auto-correlation (delay 75 -> like in original sequest algorithm)
    /// ! Uses 10 as number of windows for window normalization (like in original sequest algorithm)    
    let private spectrumToIntensityArrayMinusAutoCorrelation (lowerScanLimit,upperScanLimit) (spectrum:PeakArray<_>) =
        let si  = PeakArray.peaksToNearestUnitDaltonBinVector spectrum lowerScanLimit upperScanLimit
        let nsi = windowNormalizeIntensities si 10    
        let nsi' = autoCorrelation 75 nsi
        (nsi - nsi')  


    /// Calculates the Cross-Correlation of p_nis and ms_nis
    let private calcXCorr (p_nis:Vector<float>) (ms_nis:Vector<float>) =
        let tmp = Vector.dot p_nis ms_nis 
        if tmp < 0. then 0. else tmp 


    /// Converts the fragment ion ladders to a theoretical Sequestlike spectrum at a given charge state. 
    /// Subsequently, the spectrum is binned to the nearest mz bin (binwidth = 1 Da). Filters out peaks 
    /// that are not within the scanLimits.
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs predictOf scanlimits chargeState possiblePeptideInfos

    /// Calculates the SequestLike Scores for all peptides (possiblePeptideInfos).
    let calcSequestLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) (scanlimits) (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =                             
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // float charge
        let fCharge = float chargeState
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge
        //
        let ides = 
            possiblePeptideInfos 
            |> List.fold (fun acc lookUpResult -> 
                              let sequence = lookUpResult.BioSequence        
                              let ionSeries = calcIonSeries massfunction sequence    
                              //predicted  normalized intensity array (spectrum) 
                              let p_nis = predictOf scanlimits fCharge ionSeries 
                              let xcorr = calcXCorr p_nis ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass  (List.length sequence) xcorr nan nan
                              
                              let revPeptide_decoy = sequence |> List.rev
                              let ionSeries_decoy  = calcIonSeries massfunction revPeptide_decoy
                              let p_nis_decoy      = predictOf scanlimits fCharge ionSeries_decoy 
                              let xcorr_decoy      = calcXCorr p_nis_decoy ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass revPeptide_decoy.Length xcorr_decoy nan nan
                              targetScore::decoyScore::acc  
                     ) []
        calcNormDeltaBestToRest (ides |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext

    /// Calculates the SequestLike Scores for all theoretical spectra.    
    let calcSequestScore scanlimits (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<Vector<float>>>) spectrumID = 
        let fCharge = float chargeState
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge        

        let ides = 
            theoreticalSpectra 
            |> List.fold (fun acc theoreticalSpectrum -> 
                              let lookUpResult = theoreticalSpectrum.LookUpResult
                              let sequence = lookUpResult.BioSequence
                              let p_nis = theoreticalSpectrum.TheoSpec 
                              let xcorr = calcXCorr p_nis ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan nan
                              
                              let p_nis_decoy = theoreticalSpectrum.DecoyTheoSpec
                              let xcorr_decoy = calcXCorr p_nis_decoy  ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan nan
                              targetScore::decoyScore::acc  
                     ) []

        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext
        
    /// Calculates the sequestLike Scores for all theoretical spectra. Implemented using Async parallel. 
    let calcSequestScoreParallel scanlimits (spectrum:PeakArray<_>) scanTime chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<Vector<float>>>) spectrumID = 
        let fCharge = float chargeState
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge     
        let computeScores theoreticalSpectrum =
            let lookUpResult = theoreticalSpectrum.LookUpResult
            let sequence = lookUpResult.BioSequence
            let p_nis = theoreticalSpectrum.TheoSpec 
            let xcorr = calcXCorr p_nis ms_nis
            let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan nan
                              
            let p_nis_decoy      = theoreticalSpectrum.DecoyTheoSpec
            let xcorr_decoy      = calcXCorr p_nis_decoy ms_nis
            let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false scanTime lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan nan
            [targetScore;decoyScore]
        
        let ides = 
            [for i in theoreticalSpectra ->  async {return computeScores i}]
            |> Async.Parallel 
            |> Async.RunSynchronously
            |> List.concat    
        calcNormDeltaBestToRest (ides  |> List.sortBy (fun sls -> - sls.Score))        
        |> calcNormDeltaNext