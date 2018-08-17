namespace BioFSharp.Mz

open BioFSharp
open SearchDB


module SearchEngineResult = 
    
    type SearchEngine = 
        | AndromedaLike = 0 
        | SEQUESTLike   = 1

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
