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
    let getTheoSpec (predictTheoSpec: float*float -> float -> PeakFamily<TaggedMass.TaggedMass> list -> 'a) scanlimits chargeState ((lookUpResult,ionSeries):LookUpResult<AminoAcids.AminoAcid>*Fragmentation.FragmentMasses) =
        let theoSpec = 
            predictTheoSpec scanlimits (float chargeState) ionSeries.TargetMasses
        let theoSpecDecoy = 
            predictTheoSpec scanlimits (float chargeState) ionSeries.DecoyMasses
        createTheoreticalSpectrum lookUpResult theoSpec theoSpecDecoy
                         
    ///
    let getTheoSpecs (predictTheoSpec: float*float -> float -> PeakFamily<TaggedMass.TaggedMass> list -> 'a) scanlimits chargeState (possiblePeptideInfos:(LookUpResult<AminoAcids.AminoAcid>*Fragmentation.FragmentMasses) list) =
        possiblePeptideInfos 
        |> List.fold (fun acc lookUpResult -> 
                        getTheoSpec predictTheoSpec scanlimits chargeState lookUpResult :: acc
                    ) []

    