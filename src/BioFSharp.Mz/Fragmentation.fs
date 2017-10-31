namespace BioFSharp.Mz

module Fragmentation =

    open System
    open FSharp.Care
    open FSharp.Care.Collections
    open BioFSharp
    open BioFSharp.IO

    let waterLossSet = set [AminoAcids.Ser;AminoAcids.Thr;AminoAcids.Glu;AminoAcids.Asp;]
    let aminoLossSet = set [AminoAcids.Arg;AminoAcids.Lys;AminoAcids.Gln;AminoAcids.Asn;]

    let isWaterLoss a =
        waterLossSet.Contains(a)

    let isAminoLoss a =
        aminoLossSet.Contains(a)

    //TODO: Implement neutral loss logic 
    //let isNeutralLoss a = 
    //    match a with 
    //    | AminoAcids.Mod (aa,_) -> false
    //    | _                     -> false
    
    
    type FragmentMasses = {
        TargetMasses : PeakFamily<TaggedMass.TaggedMass> list
        DecoyMasses  : PeakFamily<TaggedMass.TaggedMass> list
        }
        
    let createFragmentMasses targetMasses decoyMasses= {
        TargetMasses = targetMasses 
        DecoyMasses  = decoyMasses
        }
    
    [<AutoOpenAttribute>]
    module private BioList =

        open BioFSharp.Mz.Peaks
        open BioFSharp.Mz.Ions
        ///
        let private calcBorYIonFragMass (massfunction:Formula.Formula -> float) acc aa = 
            acc + massfunction (AminoAcids.formula aa)
        
        ///
        let private calcAIonFragMass (massfunction:Formula.Formula -> float) acc aa =
            acc + massfunction (AminoAcids.formula aa) - (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.CO)) 
        
        ///
        let private calcCIonFragMass (massfunction:Formula.Formula -> float) acc aa =
            acc + massfunction (AminoAcids.formula aa) + (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.NH3))

        ///
        let private calcXIonFragMass (massfunction:Formula.Formula -> float) acc aa =
            acc + massfunction (AminoAcids.formula aa) + (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.CO)) 

        ///
        let private calcZIonFragMass (massfunction:Formula.Formula -> float) acc aa =
            acc + massfunction (AminoAcids.formula aa) - (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.NH3))
       
        ///
        let private createBIonTaggedMass (massfunction:Formula.Formula -> float) acc aa = 
            TaggedMass.createTaggedMass IonTypeFlag.B (calcBorYIonFragMass massfunction acc aa)

        ///
        let private createAIonTaggedMass (massfunction:Formula.Formula -> float) acc aa =
            TaggedMass.createTaggedMass IonTypeFlag.A (calcAIonFragMass massfunction acc aa)

        ///
        let private createCIonTaggedMass (massfunction:Formula.Formula -> float) acc aa =
            TaggedMass.createTaggedMass IonTypeFlag.C (calcCIonFragMass massfunction acc aa)

        ///
        let private createYIonTaggedMass (massfunction:Formula.Formula -> float) acc aa = 
            TaggedMass.createTaggedMass IonTypeFlag.Y (calcBorYIonFragMass massfunction acc aa)
        
        ///
        let private createXIonTaggedMass (massfunction:Formula.Formula -> float) acc aa =
            TaggedMass.createTaggedMass IonTypeFlag.X (calcXIonFragMass massfunction acc aa)

        ///
        let private createZIonTaggedMass (massfunction:Formula.Formula -> float) acc aa =
            TaggedMass.createTaggedMass IonTypeFlag.Z (calcZIonFragMass massfunction acc aa)

        ///
        let private peptideLadderElementOf (mainIonF: (Formula.Formula -> float) -> float -> AminoAcids.AminoAcid -> TaggedMass.TaggedMass) massfunction waterloss aminoloss acc (aa:AminoAcids.AminoAcid) = 
            //        
            let mainPeak = mainIonF massfunction acc aa    
            ///
            let lossNH3Ion  = 
                if aminoloss then 
                    let mass = (mainPeak.Mass - (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.NH3))) 
                    Some (TaggedMass.createTaggedNH3Loss mainPeak.Iontype mass)
                else  
                    None
            ///
            let lossH2OIon  = 
                if waterloss then 
                    let mass = (mainPeak.Mass - (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.H2O))) 
                    Some (TaggedMass.createTaggedH2OLoss mainPeak.Iontype mass)
                else  
                    None
            ///
            let dependentPeaks = 
                [lossNH3Ion;lossH2OIon;]//lossNeutral] 
                |> List.fold (fun acc annotPeak -> 
                                match annotPeak with 
                                | Some peak -> peak :: acc
                                | None      -> acc
                                ) []
            createPeakFamily mainPeak dependentPeaks

        ///
        let abcfragmentMassesOf (massfunction:Formula.Formula -> float) (ionSeries: Ions.IonTypeFlag) (aal:AminoAcids.AminoAcid list)  = 
            ///
            let rec series aminoList waterLoss aminoLoss fragMasses acc =
                match aminoList with
                | aa::aa'::rest    -> 
                    let waterLoss'  = waterLoss || (isWaterLoss aa)
                    let aminoLoss'  = aminoLoss || (isAminoLoss aa)
                    let bPeakFam = 
                        peptideLadderElementOf createBIonTaggedMass massfunction waterLoss' aminoLoss' acc aa
                    let aPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.A) then
                            peptideLadderElementOf createAIonTaggedMass massfunction waterLoss' aminoLoss' acc aa 
                            |> Some 
                        else 
                            None
                    let cPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.C) then
                            peptideLadderElementOf createCIonTaggedMass massfunction waterLoss' aminoLoss' acc aa 
                            |> Some 
                        else 
                            None
                    let peakFamilies =  
                        if ionSeries.HasFlag(IonTypeFlag.B) then 
                            [aPeakFam;Some bPeakFam;cPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                        else 
                            [aPeakFam;cPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                    series (aa'::rest) waterLoss' aminoLoss' (peakFamilies@fragMasses) bPeakFam.MainPeak.Mass
                | aa::rest     ->
                    let waterLoss'  = waterLoss || (isWaterLoss aa)
                    let aminoLoss'  = aminoLoss || (isAminoLoss aa)
                    let bPeakFam = 
                        peptideLadderElementOf createBIonTaggedMass massfunction waterLoss' aminoLoss' acc aa
                    let aPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.A) then
                            peptideLadderElementOf createAIonTaggedMass massfunction waterLoss' aminoLoss' acc aa 
                            |> Some 
                        else 
                            None
                    let peakFamilies =  
                        if ionSeries.HasFlag(IonTypeFlag.B) then 
                            [aPeakFam;Some bPeakFam;]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                        else 
                            [aPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                    series [] waterLoss' aminoLoss' (peakFamilies@fragMasses) bPeakFam.MainPeak.Mass
                | []          -> 
                    fragMasses
            series aal false false [] 0.0 
            |> List.rev

        ///
        let xyzfragmentMassesOf (massfunction:Formula.Formula -> float) (ionSeries: Ions.IonTypeFlag) (aal:AminoAcids.AminoAcid list)  = 
            ///
            let rec series aminoList waterLoss aminoLoss fragMasses acc =
                match aminoList with
                | aa::aa'::rest    -> 
                    let waterLoss'  = waterLoss || (isWaterLoss aa)
                    let aminoLoss'  = aminoLoss || (isAminoLoss aa)
                    let yPeakFam = 
                        peptideLadderElementOf createYIonTaggedMass massfunction waterLoss' aminoLoss' acc aa
                    let xPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.X) then
                            peptideLadderElementOf createXIonTaggedMass massfunction waterLoss' aminoLoss' acc aa |> Some 
                        else 
                            None
                    let zPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.Z) then
                            peptideLadderElementOf createZIonTaggedMass massfunction waterLoss' aminoLoss' acc aa |> Some 
                        else 
                            None
                    let peakFamilies =  
                        if ionSeries.HasFlag(IonTypeFlag.Y) then 
                            [zPeakFam;Some yPeakFam;xPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                        else 
                            [zPeakFam;xPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                    series (aa'::rest) waterLoss' aminoLoss' (peakFamilies@fragMasses) yPeakFam.MainPeak.Mass
                | aa::rest     ->
                    let waterLoss'  = waterLoss || (isWaterLoss aa)
                    let aminoLoss'  = aminoLoss || (isAminoLoss aa)
                    let yPeakFam = 
                        peptideLadderElementOf createYIonTaggedMass massfunction waterLoss' aminoLoss' acc aa
                    let zPeakFam = 
                        if ionSeries.HasFlag(IonTypeFlag.Z) then
                            peptideLadderElementOf createZIonTaggedMass massfunction waterLoss' aminoLoss' acc aa |> Some 
                        else 
                            None
                    let peakFamilies =  
                        if ionSeries.HasFlag(IonTypeFlag.Y) then 
                            [zPeakFam;Some yPeakFam;]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                        else 
                            [zPeakFam]
                            |> List.fold (fun acc peakFam -> 
                                    match peakFam with 
                                    | Some peakF -> peakF :: acc
                                    | None      ->  acc
                                    ) []
                    series [] waterLoss' aminoLoss' (peakFamilies@fragMasses) yPeakFam.MainPeak.Mass
                | []          -> 
                    fragMasses
            let ySeries = (aal |> List.rev)
            series ySeries false false [] (massfunction Formula.Table.H2O) 

    module Series = 
                    
        //type NTerminalSeries = (Formula.Formula -> float) -> AminoAcids.AminoAcid list -> PeakFamily<TaggedMass.TaggedMass> list
        //type CTerminalSeries = (Formula.Formula -> float) -> AminoAcids.AminoAcid list -> PeakFamily<TaggedMass.TaggedMass> list 
        ///// 
        //let abcOfBioList :NTerminalSeries = 
        //    fun massfunction aal -> BioList.abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.B + Ions.IonTypeFlag.C) aal 
        
        /// 
        let abcOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            BioList.abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.B + Ions.IonTypeFlag.C) aal 

        ///
        let abOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) =  
            (abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.B) aal)

        ///
        let acOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.C) aal

        ///
        let bcOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.B + Ions.IonTypeFlag.C) aal

        ///
        let aOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A) aal

        ///
        let bOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.B) aal

        ///
        let cOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.C) aal

        ///
        let xyzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Y + Ions.IonTypeFlag.Z) aal

        ///
        let xyOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Y) aal

        ///
        let xzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Z) aal

        ///
        let yzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Y + Ions.IonTypeFlag.Z) aal

        ///
        let xOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X) aal

        ///
        let yOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Y) aal

        ///
        let zOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Z) aal

        let inline fragmentMasses (nTerminalSeries) (cTerminalSeries) (massFunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            let targetMasses = 
                let nTerm = nTerminalSeries massFunction aal
                let cTerm = cTerminalSeries massFunction aal
                nTerm@cTerm
            let decoyMasses = 
                let nTerm = nTerminalSeries massFunction (aal |> List.rev)
                let cTerm = cTerminalSeries massFunction (aal |> List.rev)
                nTerm@cTerm
            createFragmentMasses targetMasses decoyMasses        
 
 //    let imoniumIons (rawMass:List<float>) (label : Formula.Formula -> Formula.Formula) = 
    //        let currentCO = massDiffAX_CO label 
    //        rawMass |> List.map (fun n ->  n - currentCO)


    //TODO: Integrate in fragmentation method 
    //let lossNeutral = 
    //    if BioFSharp.Mz.Fragmentation.isNeutralLoss f then 
    //        match f with 
    //        | AminoAcids.Mod(aa, modiL) ->  
    //        //TODO: Implement neutral loss logic    
    //            None
    //        //    let mass = (acc + currentMass - (massfunction (AminoAcids.isotopicLabelFunc f (modiL.Head.Modify Formula.emptyFormula)) ) )
    //        //    Some (createTag (true,(NeutralLoss ionSeries)) mass)
    //        | _ -> None
    //    else None
