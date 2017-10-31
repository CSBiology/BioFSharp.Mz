namespace BioFSharp.Mz

open BioFSharp
open SearchDB
open Fragmentation
open TheoreticalSpectra
open FSharp.Care
open SearchEngineResult

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra.Double

module SequestLike =

    /// normalize the intensities within a window to maximum of the window
    /// Attention shortens the array  (cuts)
    let windowNormalizeIntensities (intensities:LinearAlgebra.Vector<float>) (numberOfWindows:int) =
        // finds max within range
        let rec findMax (array:LinearAlgebra.Vector<float>) (cMax:float) (lowerLimit:int) (counter:int)  =
            if counter < lowerLimit then
                cMax
            else
                let nMax = max array.[counter] cMax
                findMax array nMax lowerLimit (counter - 1)
        // sqrt intensities and normalize to maximum within lower-upper limit
        let normToMaxSqrtInPlace (array:LinearAlgebra.Vector<float>) (lowerLimit:int) (upperLimit:int)  =
            let nMax = sqrt (findMax array 0.0 lowerLimit upperLimit)
            for i = lowerLimit to upperLimit do
                if nMax > 0. then                
                    array.[i] <- (sqrt array.[i]) / nMax 
                else
                    array.[i] <- 0.
        let windowSize =  (intensities.Count / numberOfWindows)    
        let tmpIntensities =            
            DenseVector.Create( (windowSize*numberOfWindows) ,(fun i -> intensities.[i]))                    
        for i = 1 to numberOfWindows do
            //printfn "window: %i lower: %i counter: %i " i (windowSize * (i - 1)) (windowSize * i - 1)
            normToMaxSqrtInPlace tmpIntensities (windowSize * (i - 1)) (windowSize * i - 1) 
        tmpIntensities
    
    ///
    let predictIntensitySimpleModel (ionType:Ions.IonTypeFlag) (charge:float) =
        match ionType with 
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Unknown     -> 1.  / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Diagnostic  -> 1.  / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Neutral     -> 0.6 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Immonium    -> 0.6 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.lossNH3     -> 0.2 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.lossH2O     -> 0.2 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.A           -> 0.2 / charge   
        | iontype when ionType.HasFlag Ions.IonTypeFlag.B           -> 1.  / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.C           -> 0.2 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.X           -> 0.2 / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Y           -> 1.  / charge
        | iontype when ionType.HasFlag Ions.IonTypeFlag.Z           -> 0.2 / charge
        | _                                                         -> 0.2 / charge

    ///
    let predictOf (maxcharge:float) (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
        let predictPeak charge (taggedMass: TaggedMass.TaggedMass) = 
            Peak((Mass.toMZ taggedMass.Mass charge), (predictIntensitySimpleModel (taggedMass.Iontype) charge))
        let rec recloop (ions) (charge) (fragments:PeakFamily<TaggedMass.TaggedMass> list) =
            match fragments with
            | fragment::rest ->  
                let mainPeak = 
                    predictPeak charge fragment.MainPeak
                let dependendPeaks =
                    fragment.DependentPeaks
                    |> List.map (predictPeak charge)
                recloop (mainPeak::dependendPeaks@ions) charge rest
            | []            -> ions
            
        [| for z = 1. to maxcharge do
                yield! recloop [] z fragments |]

    ///
    let shiftedVectorSum (plusMinusMaxDelay:int) (vector:DenseVector) =
        let shifted (vector:DenseVector) (tau:int) =
            vector
            |> LinearAlgebra.Vector.mapi
                (fun i x ->
                    let index = i - tau
                    if (index < 0) || (index > vector.Count - 1) then 
                        0.
                    else
                        vector.[index] )
        let rec accumVector (accum) (state:int) (max:int) =
            if state = max then
                accum
            else
                accumVector (accum + (shifted vector state)) (state - 1) (max)
        let emtyVector = DenseVector (Array.zeroCreate vector.Count)
        let plus  = accumVector emtyVector (plusMinusMaxDelay) (1)
        let minus = accumVector emtyVector (-1) (-plusMinusMaxDelay)
        (plus + minus)
        |> LinearAlgebra.Vector.map (fun x ->  x / (float plusMinusMaxDelay * 2.))
       
    ///
    let createShiftedMatrix (plusMinusMaxDelay:int) (array:float[]) =
        let colNumber = array.Length
        Array2D.init (plusMinusMaxDelay * 2 + 1) colNumber
            (fun i ii ->
                let ni = (i - plusMinusMaxDelay)
                let index = ii - ni
                if (index < 0) || (index > colNumber - 1) then 
                    0.
                else
                    array.[index] )

    /// Amino acid sequence (peptide) to sequest-like predicted intensity array
    let peaksToNormalizedIntensityArray (lowerScanLimit,upperScanLimit) charge ionSeries =         
        let psi  = predictOf charge ionSeries  
        let npsi = PeakArray.peaksToNearestUnitDaltonBinVector psi lowerScanLimit upperScanLimit        
        npsi

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
    let spectrumToIntensityArrayMinusAutoCorrelation (lowerScanLimit,upperScanLimit) (spectrum:PeakArray<_>) =
        let si  = PeakArray.peaksToNearestUnitDaltonBinVector spectrum lowerScanLimit upperScanLimit
        let nsi = windowNormalizeIntensities si 10    
        let nsi' = shiftedVectorSum 75 nsi
        (nsi - nsi')  


    /// Calculates sequest-like deltaCN score
    ///  (Xcorr(top hit) - Xcorr(n)) ÷ Xcorr(top hit). Thus, the deltaCn for the top hit is
    ///  (Xcorr(top hit) - Xcorr(top hit)) ÷ Xcorr(top hit) = 0.
    let calcDeltaCN (sourceList:SearchEngineResult<'a> list) =
        match sourceList with
        | h1::rest -> 
            sourceList
            |> List.map
                (fun sls ->
                    let deltaCN = (h1.Score - sls.Score) / h1.Score
                    { sls with DeltaCN = deltaCN } )      
        | []       -> []

    ///
    let calcSequestLikeScoresRevDecoy calcIonSeries (massfunction:Formula.Formula -> float) (scanlimits) (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) spectrumID =                             
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
                              let p_nis = peaksToNormalizedIntensityArray scanlimits fCharge ionSeries 
                              let xcorr = p_nis * ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass  (List.length sequence) xcorr nan
                              
                              let revPeptide_decoy = sequence |> List.rev
                              let ionSeries_decoy  = calcIonSeries massfunction revPeptide_decoy
                              let p_nis_decoy      = peaksToNormalizedIntensityArray scanlimits fCharge ionSeries_decoy 
                              let xcorr_decoy      = p_nis_decoy * ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass revPeptide_decoy.Length xcorr_decoy nan 
                              targetScore::decoyScore::acc  
                     ) []
        calcDeltaCN (ides |> List.sortBy (fun sls -> - sls.Score))        

    ///
    let getTheoSpecs scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*FragmentMasses>) =
        TheoreticalSpectra.getTheoSpecs peaksToNormalizedIntensityArray scanlimits chargeState possiblePeptideInfos

    ///          
    let calcSequestScore scanlimits (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<LinearAlgebra.Vector<float>>>) spectrumID = 
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
                              let xcorr = p_nis * ms_nis
                              let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan
                              
                              let p_nis_decoy      = theoreticalSpectrum.DecoyTheoSpec
                              let xcorr_decoy      = p_nis_decoy * ms_nis
                              let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan 
                              targetScore::decoyScore::acc  
                     ) []

        calcDeltaCN (ides  |> List.sortBy (fun sls -> - sls.Score))        
            
    ///          
    let calcSequestScoreParallel scanlimits (spectrum:PeakArray<_>) chargeState isolationWindowTargetMz (theoreticalSpectra:list<TheoreticalSpectrum<LinearAlgebra.Vector<float>>>) spectrumID = 
        let fCharge = float chargeState
        // measured normailzed intensity array (spectrum) minus auto-correlation
        let ms_nis =  spectrumToIntensityArrayMinusAutoCorrelation scanlimits spectrum
        // measured mass
        let ms_mass = Mass.ofMZ isolationWindowTargetMz fCharge     
        let computeScores theoreticalSpectrum =
            let lookUpResult = theoreticalSpectrum.LookUpResult
            let sequence = lookUpResult.BioSequence
            let p_nis = theoreticalSpectrum.TheoSpec 
            let xcorr = p_nis * ms_nis
            let targetScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod true lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr nan
                              
            let p_nis_decoy      = theoreticalSpectrum.DecoyTheoSpec
            let xcorr_decoy      = p_nis_decoy * ms_nis
            let decoyScore = createSearchEngineResult SearchEngineResult.SearchEngine.SEQUESTLike spectrumID lookUpResult.ModSequenceID lookUpResult.PepSequenceID lookUpResult.GlobalMod false lookUpResult.StringSequence chargeState isolationWindowTargetMz ms_mass lookUpResult.Mass sequence.Length xcorr_decoy nan 
            [targetScore;decoyScore]
        
        let ides = 
            [for i in theoreticalSpectra ->  async {return computeScores i}]
            |> Async.Parallel 
            |> Async.RunSynchronously
            |> List.concat    
        calcDeltaCN (ides  |> List.sortBy (fun sls -> - sls.Score))        
    