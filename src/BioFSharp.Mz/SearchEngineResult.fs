namespace BioFSharp.Mz

open BioFSharp
open SearchDB


module SearchEngineResult = 
    
    type SearchEngine = 
        | AndromedaLike 
        | SEQUESTLike   
        | XTandemLike

    type SearchEngineResult<'a> = {
        SearchEngine          : SearchEngine
        SpectrumID            : string;     
        ModSequenceID         : int
        PepSequenceID         : int;
        GlobalMod             : int;
        IsTarget              : bool;
        ScanTime              : float; 
        StringSequence        : string;
        PrecursorCharge       : int;
        PrecursorMZ           : float
        MeasuredMass          : float;
        TheoMass              : float;            
        PeptideLength         : int;
        Score                 : 'a;
        NormDeltaBestToRest   : float;
        NormDeltaNext         : float;}
    
    ///
    let createSearchEngineResult searchEngine spectrumID modSequenceID pepSequenceID globalMod isTarget scanTime stringSequence precursorCharge precursorMZ measuredMass theoMass peptideLength xcorr normDeltaBestToRest normDeltaNext = 
        { SearchEngine=searchEngine;SpectrumID = spectrumID; ModSequenceID = modSequenceID; PepSequenceID = pepSequenceID;ScanTime = scanTime ; GlobalMod = globalMod; IsTarget = isTarget; StringSequence = stringSequence; 
          PrecursorCharge = precursorCharge; PrecursorMZ=precursorMZ;  MeasuredMass = measuredMass; TheoMass = theoMass; PeptideLength = peptideLength; Score = xcorr; NormDeltaBestToRest = normDeltaBestToRest;NormDeltaNext = normDeltaNext }

    /// Calculates sequest-like delta normalized by best score to the best score.
    ///  (Xcorr(top hit) - Xcorr(n)) / Xcorr(top hit). Thus, the deltaCn for the top hit is
    ///  (Xcorr(top hit) - Xcorr(top hit)) / Xcorr(top hit) = 0.
    /// if the best Score equals 0. this function returns returns 1 for every PSM
    let calcNormDeltaBestToRest (sourceList:SearchEngineResult<'a> list) =
        match sourceList with
        | h1::rest -> 
            if h1.Score <= 0. then sourceList |> List.map (fun sls -> {sls with NormDeltaBestToRest = 1.})
            else
            sourceList
            |> List.map
                (fun sls ->
                    let normDeltaBestToRest  = (h1.Score - sls.Score) / h1.Score
                    { sls with NormDeltaBestToRest = normDeltaBestToRest } )      
        | []       -> []

    // Iterates over the score ranked PSMs and computes the score difference between adjacent
    // PSMs normalized by the Score of the best ranked PSM.
    /// if the best Score equals 0., this function returns returns 0 for every PSM.
    let calcNormDeltaNext (sourceList:SearchEngineResult<'a> list) =        
        let rec loop normF acc l = 
            match l with 
            | hLast::[] ->
                {hLast with NormDeltaNext = 0.}::acc
                |> List.rev 
            | hi::hii -> 
                let normDeltaNext = (hi.Score - hii.[0].Score) / normF 
                loop normF ({hi with NormDeltaNext = normDeltaNext}::acc) hii 
        match sourceList with
        | h1::rest -> 
            let normFactor = h1.Score
            if normFactor <= 0. then sourceList |> List.map (fun sls -> {sls with NormDeltaNext = 0.})
            else
            loop normFactor [] sourceList                  
        | []       -> []
