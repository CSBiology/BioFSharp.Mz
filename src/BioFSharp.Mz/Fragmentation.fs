namespace BioFSharp.Mz

module Fragmentation =

    open System
    open FSharpAux
    open BioFSharp
    open BioFSharp.IO

    /// Set of amino acids that are prone to cleave off an H2O molecule of their side chain uppon fragmentation.
    let waterLossSet = set [AminoAcids.Ser;AminoAcids.Thr;AminoAcids.Glu;AminoAcids.Asp;]

    /// Set of amino acids that are prone to cleave off an NH2 molecule of their side chain uppon fragmentation.
    let aminoLossSet = set [AminoAcids.Arg;AminoAcids.Lys;AminoAcids.Gln;AminoAcids.Asn;]

    /// Returns true if the amino acid "a" is prone to cleave off an H2O molecule of their side chain uppon fragmentation.
    let isWaterLoss a =
        waterLossSet.Contains(a)

    /// Returns true if the amino acid "a" is prone to cleave off an NH2 molecule of their side chain uppon fragmentation.
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
    
    [<AutoOpen>]
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


        /// Returns the masses of a, b and c series, specified by the ionSeries parameter. The mass accuracy is determined by the massfunction applied.        
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

        /// Returns the masses of x, y and z series, specified by the ionSeries parameter. The mass accuracy is determined by the massfunction applied.
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
                    
        /// Retunrs the a b and c series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let abcOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.B + Ions.IonTypeFlag.C) aal 

        /// Returns the a and b series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let abOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) =  
            (abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.B) aal)

        /// Returns the a and c series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let acOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A + Ions.IonTypeFlag.C) aal

        /// Returns the b and c series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let bcOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.B + Ions.IonTypeFlag.C) aal

        /// Returns the a series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let aOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.A) aal

        /// Returns the b series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let bOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.B) aal

        /// Returns the c series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let cOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            abcfragmentMassesOf massfunction (Ions.IonTypeFlag.C) aal

        /// Returns the x y and z series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let xyzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Y + Ions.IonTypeFlag.Z) aal

        /// Returns the x and y series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let xyOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Y) aal

        /// Returns the x and z series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let xzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X + Ions.IonTypeFlag.Z) aal

        /// Returns the y and z series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let yzOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Y + Ions.IonTypeFlag.Z) aal

        /// Returns the x series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let xOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.X) aal

        /// Returns the y series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let yOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Y) aal

        /// Returns the z series of the given amino acids list. The mass accuracy is determined by the massfunction applied.
        let zOfBioList (massfunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
            xyzfragmentMassesOf massfunction (Ions.IonTypeFlag.Z) aal

        /// Returns the fragment masses of the amino acid sequence specified by aal. The ionseries are specified by functions
        /// "nTerminalSeries" and "cTerminalSeries". The mass accuracy is determined by the massfunction applied.
        let inline fragmentMasses nTerminalSeries cTerminalSeries (massFunction:Formula.Formula -> float) (aal:AminoAcids.AminoAcid list) = 
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
