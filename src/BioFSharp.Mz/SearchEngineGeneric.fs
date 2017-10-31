namespace BioFSharp.Mz

open BioFSharp
open SearchDB
open TheoreticalSpectra

module SearchEngineGeneric = 

    module OrderedCache =

        /////                
        let private getTheoSpecsWithMem predictTheoSpec (spectrumCache: Cache.Cache<int64,TheoreticalSpectrum<'a> list>) scanlimits chargeState (possiblePeptideInfos:list<LookUpResult<AminoAcids.AminoAcid>*Fragmentation.FragmentMasses>) =
            if possiblePeptideInfos.IsEmpty then []
            else 
            let fCharge = chargeState |> float
            let lowerMassAsInt = possiblePeptideInfos |> List.minBy (fun (lookUpRes,frag) -> lookUpRes.RoundedMass) |> fun (lookUpRes,frag) -> lookUpRes.RoundedMass
            let upperMassAsInt = possiblePeptideInfos |> List.maxBy (fun (lookUpRes,frag) -> lookUpRes.RoundedMass) |> fun (lookUpRes,frag) -> lookUpRes.RoundedMass
            match Cache.containsItemsBetween spectrumCache (lowerMassAsInt,upperMassAsInt) with
            // The cache contains items that lie within the search range
            | Some (lowerMassIdx,upperMassIdx)  -> 
                let upperMassInCache = spectrumCache.Keys.[upperMassIdx]
                /// generate the spectra missing in the cache 
                possiblePeptideInfos
                |> List.fold (fun acc (lookUpRes,frag) -> 
                                match lookUpRes.RoundedMass >= lowerMassAsInt && lookUpRes.RoundedMass <= upperMassInCache with
                                | true  -> 
                                    acc
                                | false ->
                                    let theoSpec = 
                                        predictTheoSpec scanlimits (float chargeState) frag.TargetMasses
                                    let decoySeq = lookUpRes.BioSequence |> List.rev
                                    let theoSpecDecoy = 
                                        predictTheoSpec scanlimits (float chargeState) frag.DecoyMasses
                                    (createTheoreticalSpectrum lookUpRes theoSpec theoSpecDecoy)::acc
                                ) [] 
                /// add the missing spectra to the cache                
                |> List.groupBy (fun theoSpec -> theoSpec.LookUpResult.RoundedMass)
                |> List.map (fun x -> Cache.addItem spectrumCache x)
                |> ignore                 
                /// look up theoretical spectra in the updated cache 
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
            // The cache contains no items of this mass that lie within the search range.
            | None ->
                /// generate the spectra missing in the cache 
                let missingSpectra = 
                    possiblePeptideInfos
                    |> List.fold (fun acc (lookUpRes,frag) -> 
                                        let theoSpec = 
                                            predictTheoSpec scanlimits (float chargeState) frag.TargetMasses
                                        let decoySeq = lookUpRes.BioSequence |> List.rev
                                        let theoSpecDecoy = 
                                            predictTheoSpec scanlimits (float chargeState) frag.DecoyMasses
                                        (createTheoreticalSpectrum lookUpRes theoSpec theoSpecDecoy)::acc
                                    ) [] 
                /// add the missing spectra to the cache
                missingSpectra 
                |> List.groupBy (fun theoSpec -> theoSpec.LookUpResult.RoundedMass)
                |> List.map (fun x -> Cache.addItem spectrumCache x)          
                |> ignore 
                // return spectra added to the cache.
                missingSpectra   
                           


        /// 
        let generateTheoSpectra calcIonSeries (massfunction:Formula.Formula -> float) (lookUpF: float -> float -> LookUpResult<AminoAcids.AminoAcid> list) (lookUpCache: Cache.Cache<int64,((LookUpResult<AminoAcids.AminoAcid>*Fragmentation.FragmentMasses) list)>) 
                (andromedaCache: Cache.Cache<int64,TheoreticalSpectrum<'a> list>) (sequestCache: Cache.Cache<int64,TheoreticalSpectrum<'b> list>)  chargeState scanlimits
                    (maxMemory:int64)  lowerMass upperMass = 
            ///
            let memoryUsage = System.GC.GetTotalMemory(false)
            if memoryUsage > maxMemory then
                lookUpCache.Clear()
                andromedaCache.Clear()
                sequestCache.Clear()
            ///
            let lookUpResults = SearchDB.getPeptideLookUpWithMemBy calcIonSeries massfunction lookUpF lookUpCache lowerMass upperMass
            let andromedaResults = getTheoSpecsWithMem AndromedaLike.predictOf andromedaCache scanlimits chargeState lookUpResults
            let sequestResults   = getTheoSpecsWithMem SequestLike.peaksToNormalizedIntensityArray sequestCache scanlimits chargeState lookUpResults
            andromedaResults,sequestResults