namespace BioFSharp.Mz

open BioFSharp
open SearchDB


module SearchEngineResult = 
    
    type SearchEngine = 
        | AndromedaLike = 0 
        | SEQUESTLike   = 1

    type SearchEngineResult<'a> = {
        SearchEngine    : SearchEngine
        SpectrumID      : string;     
        ModSequenceID   : int
        PepSequenceID   : int;
        GlobalMod       : int;
        IsTarget        : bool;
        StringSequence  : string;
        PrecursorCharge : int;
        PrecursorMZ     : float
        MeasuredMass    : float;
        TheoMass        : float;            
        PeptideLength   : int;
        Score           : 'a;
        DeltaCN         : float; }
    
    ///
    let createSearchEngineResult searchEngine spectrumID modSequenceID pepSequenceID  globalMod isTarget peptide precursorCharge precursorMZ measuredMass theoMass peptideLength xcorr deltaCN = 
        { SearchEngine=searchEngine;SpectrumID = spectrumID; ModSequenceID = modSequenceID; PepSequenceID = pepSequenceID; GlobalMod = globalMod; IsTarget = isTarget; StringSequence = peptide; 
          PrecursorCharge = precursorCharge; PrecursorMZ=precursorMZ;  MeasuredMass = measuredMass; TheoMass = theoMass; PeptideLength = peptideLength; Score = xcorr; DeltaCN = deltaCN; }
