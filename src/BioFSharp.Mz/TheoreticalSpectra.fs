namespace BioFSharp.Mz

open BioFSharp
open SearchDB


module TheoreticalSpectra = 
    
    type TheoreticalSpectrum<'a> = {
        LookUpResult  : LookUpResult<AminoAcids.AminoAcid> 
        TheoSpec      : 'a
        DecoyTheoSpec : 'a
        }

    ///
    let createTheoreticalSpectrum lookUpResult theoSpec decoyTheoSpec = {
        LookUpResult  = lookUpResult
        TheoSpec      = theoSpec
        DecoyTheoSpec = decoyTheoSpec
        } 

    ///
    let getTheoSpec (predictTheoSpec: (Formula.Formula -> float) -> float*float -> float -> AminoAcids.AminoAcid list -> 'a) (massfunction:Formula.Formula -> float) scanlimits chargeState (lookUpResult:LookUpResult<AminoAcids.AminoAcid>) =
        let sequence = lookUpResult.BioSequence
        let theoSpec = 
            predictTheoSpec massfunction scanlimits (float chargeState) sequence
        let decoySeq = sequence |> List.rev
        let theoSpecDecoy = 
            predictTheoSpec massfunction scanlimits (float chargeState) decoySeq
        createTheoreticalSpectrum lookUpResult theoSpec theoSpecDecoy
                         
    ///
    let getTheoSpecs (predictTheoSpec: (Formula.Formula -> float) -> float*float -> float -> AminoAcids.AminoAcid list -> 'a) (massfunction:Formula.Formula -> float) scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) =
        possiblePeptideInfos 
        |> List.fold (fun acc lookUpResult -> 
                        getTheoSpec predictTheoSpec massfunction scanlimits chargeState lookUpResult :: acc
                    ) []

    