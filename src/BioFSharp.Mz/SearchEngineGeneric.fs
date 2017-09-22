namespace BioFSharp.Mz

open BioFSharp
open SearchDB
open AndromedaLike
open SequestLike
open TheoreticalSpectra
//open Theor


module SearchEngineGeneric = 

    module OrderedCache =

        /////                
        let private getTheoSpecsWithMem (massfunction:Formula.Formula -> float) (predictTheoSpec: (Formula.Formula -> float) -> float*float -> float -> AminoAcids.AminoAcid list -> 'a) (spectrumCache: Cache.Cache<int64,TheoreticalSpectrum<'a> list>) scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>>) =
            let fCharge = chargeState |> float
            let lowerMassAsInt = possiblePeptideInfos |> List.minBy (fun x -> x.RoundedMass) |> fun x -> x.RoundedMass
            let upperMassAsInt = possiblePeptideInfos |> List.maxBy (fun x -> x.RoundedMass) |> fun x -> x.RoundedMass
            match Cache.containsItemsBetween spectrumCache (lowerMassAsInt,upperMassAsInt) with
            // The cache contains items that lie within the search range
            | true, (lowerMassIdx,upperMassIdx)  -> 
                let upperMassInCache = spectrumCache.Keys.[upperMassIdx]
                let updateCache =
                    /// generate the spectra missing in the cache 
                    let missingSpectra = 
                        possiblePeptideInfos
                        |> List.fold (fun acc lookUpRes -> 
                                        match lookUpRes.RoundedMass >= lowerMassAsInt && lookUpRes.RoundedMass <= upperMassInCache with
                                        | true  -> 
                                            acc
                                        | false ->
                                            let theoSpec = 
                                                predictTheoSpec massfunction scanlimits (float chargeState) lookUpRes.BioSequence
                                            let decoySeq = lookUpRes.BioSequence |> List.rev
                                            let theoSpecDecoy = 
                                                predictTheoSpec massfunction scanlimits (float chargeState) decoySeq
                                            (createTheoreticalSpectrum lookUpRes theoSpec theoSpecDecoy)::acc
                                        ) [] 
                    /// add the missing spectra to the cache
                    let updateCache = 
                        missingSpectra 
                        |> List.groupBy (fun theoSpec -> theoSpec.LookUpResult.RoundedMass)
                        |> List.map (fun x -> Cache.addItem spectrumCache x)
                    missingSpectra
                /// look up theoretical spectra in the updated cache 
                let lookUpResults = 
                    /// assure that the lower key lies within the search range
                    if spectrumCache.Keys.[lowerMassIdx] < lowerMassAsInt then 
                        let lowerMassIdx' = lowerMassIdx + 1
                        let upperMassIdx'  = Cache.binarySearch Cache.Border.Upper spectrumCache upperMassAsInt 
                        Cache.getValuesByIdx spectrumCache (lowerMassIdx',upperMassIdx')//(lookUpCache.Count-1))
                        |> List.concat
                    else
                        let upperMassIdx'  = Cache.binarySearch Cache.Border.Upper spectrumCache upperMassAsInt 
                        //printfn "upperMassIdx = %i, specChachelength =%i " upperMassIdx' (spectrumCache.Count-1)
                        Cache.getValuesByIdx spectrumCache (lowerMassIdx,upperMassIdx')//(lookUpCache.Count-1))
                        |> List.concat
                lookUpResults
            // The cache contains no items of this mass that lie within the search range.
            | false, (lowerMassIdx,upperMassIdx) ->
                /// generate the spectra missing in the cache 
                let missingSpectra = 
                    possiblePeptideInfos
                    |> List.fold (fun acc lookUpRes -> 
                                        let theoSpec = 
                                            predictTheoSpec massfunction scanlimits (float chargeState) lookUpRes.BioSequence
                                        let decoySeq = lookUpRes.BioSequence |> List.rev
                                        let theoSpecDecoy = 
                                            predictTheoSpec massfunction scanlimits (float chargeState) decoySeq
                                        (createTheoreticalSpectrum lookUpRes theoSpec theoSpecDecoy)::acc
                                    ) [] 
                /// add the missing spectra to the cache
                let updateCache = 
                    missingSpectra 
                    |> List.groupBy (fun theoSpec -> theoSpec.LookUpResult.RoundedMass)
                    |> List.map (fun x -> Cache.addItem spectrumCache x)          
                // return spectra added to the cache.
                missingSpectra   
               
            
        ///// 
        //let generateTheoSpectra (massfunction:Formula.Formula -> float) (lookUpF: float -> float -> LookUpResult<AminoAcids.AminoAcid> list) (lookUpCache: Cache.Cache<int64,(LookUpResult<AminoAcids.AminoAcid> list)>) 
        //        (andromedaCache: Cache.Cache<int64,TheoreticalSpectrum<'a> list>) (sequestCache: Cache.Cache<int64,TheoreticalSpectrum<'b> list>)  chargeState scanlimits
        //            (maxMemory:int64)  lowerMass upperMass = 
        //    ///
        //    let memoryUsage = System.GC.GetTotalMemory(false)
        //    if memoryUsage > maxMemory then
        //        lookUpCache.Clear()
        //        andromedaCache.Clear()
        //        sequestCache.Clear()
        //    ///
        //    let lookUpResults = SearchDB.getPeptideLookUpWithMemBy lookUpF lookUpCache lowerMass upperMass
        //    let andromedaResults = getTheoSpecsWithMem massfunction AndromedaLike.predictOf andromedaCache scanlimits chargeState lookUpResults
        //    let sequestResults   = getTheoSpecsWithMem massfunction SequestLike.peptideToNormalizedIntensityArray sequestCache scanlimits chargeState lookUpResults
        //    andromedaResults,sequestResults