namespace BioFSharp.Mz


open FSharpAux
open BioFSharp.PeptideClassification
open System.Collections.Generic


module ProteinInference = 

    type internal StructuralEqualityComparer<'a when 'a : equality >() = 
        interface IEqualityComparer<'a> with
            member x.Equals(a,b) = a = b
            member x.GetHashCode a = hash a 
            
    type ProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs: string []
            PeptideSequence: 'sequence
            Class: PeptideEvidenceClass
        }

    type InferredProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs: string
            PeptideSequence: 'sequence []
            Class: PeptideEvidenceClass
        }    

    let createInferredProteinClassItem proteinIDs evidenceClass peptideSequences = 
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence = peptideSequences
            Class = evidenceClass
        }   

    let createProteinClassItem proteinIDs evidenceClass peptideSequence : ProteinClassItem<'sequence> = 
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence = peptideSequence
            Class = evidenceClass
        }
    
    type IntegrationStrictness = 
        | KeepOverlap
        | KeepSmallerGroup

    let private increaseClass c =
        match c with
        | PeptideEvidenceClass.C1a -> PeptideEvidenceClass.C1b
        | PeptideEvidenceClass.C1b -> PeptideEvidenceClass.C2a
        | PeptideEvidenceClass.C2a -> PeptideEvidenceClass.C2b
        | PeptideEvidenceClass.C2b -> PeptideEvidenceClass.C3a
        | PeptideEvidenceClass.C3a -> PeptideEvidenceClass.C3b

    let private mapClassName c =
        match c with
        | PeptideEvidenceClass.C1a -> "Class1A"
        | PeptideEvidenceClass.C1b -> "Class1B"
        | PeptideEvidenceClass.C2a -> "Class2A"
        | PeptideEvidenceClass.C2b -> "Class2B"
        | PeptideEvidenceClass.C3a -> "Class3A"
        | PeptideEvidenceClass.C3b -> "Class3B"

    /// Contains all different classes with their associated proteins

    let private proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup


    type ClassMap =
        Map<PeptideEvidenceClass,BidirectionalDictionary<string,string>>

    module private ClassMap =

        ///Creates a new empty ClassMap
        let empty() : ClassMap =
            Map.empty
            |> Map.add PeptideEvidenceClass.C1a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C1b (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C2a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C2b (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C3a (BidirectionalDictionary<string,string>())
            |> Map.add PeptideEvidenceClass.C3b (BidirectionalDictionary<string,string>())

        ///Searches the Memory Map from highest ranking to lowest ranking class for any of the given proteinIds
        let searchProts (m:ClassMap) (prots:string[]) =
            match Map.tryFindKey (fun c (d:BidirectionalDictionary<string,string>) ->
                    Array.exists (fun p -> d.ContainsKey p) prots
                ) m with
            | Some c -> c
            | None -> PeptideEvidenceClass.Unknown

        ///Adds the protein group to the given class in the memory map
        let addGroup (m:ClassMap) (c:PeptideEvidenceClass) prots =
            let name = proteinGroupToString prots
            Seq.iter (fun p -> m.[c].Add p name) prots
            m

    ///In a given class finds the overlap between the given and already added proteins and integrates them
    let integrate (m:ClassMap) c proteinGroups : ClassMap =
        let int (dict:BidirectionalDictionary<string,string>) =
            let dict = dict
            let oldGroups =
                Array.choose (dict.TryGetByKey >> (Option.map Array.ofSeq)) proteinGroups
                |> Array.concat
                |> Array.distinct
                |> Array.choose (fun pg -> Option.map (fun seqs -> seqs,pg) (dict.TryGetByValue pg))
            let newSet = Set.ofSeq proteinGroups
            let newGroup = proteinGroupToString proteinGroups
            oldGroups
            |> Array.iter (fun (x,y) -> printfn "%s" ((Seq.toArray >> proteinGroupToString) x) )
            oldGroups
            |> Array.iter (fun (seqs,group) -> 
                let oldSet = Set.ofSeq seqs
                if Set.isProperSubset newSet oldSet then 
                    dict.RemoveValue group
                    Seq.iter (fun k -> dict.Add k newGroup) proteinGroups
            ) 
            dict
        Map.map (fun ca dict ->
            if ca = c then int dict
            else dict) m
    
    let groupIntoClasses (integrationFunction: ClassMap -> PeptideEvidenceClass -> string [] -> ClassMap) (proteinClassItems: ProteinClassItem<'sequence> list ) =
        let rec loop (m:ClassMap) proteinClassItems remaining lookingForClass =
            let findAndIntegrate (m:ClassMap) (c:PeptideEvidenceClass) (proteins:string []) =
                match ClassMap.searchProts m proteins with
                | PeptideEvidenceClass.Unknown ->
                    printfn "Search did not return anything"
                    ClassMap.addGroup m c proteins
                | foundC ->
                    if foundC < c then
                        m
                    else
                        integrationFunction m foundC proteins
            match proteinClassItems with
            | (item:ProteinClassItem<'sequence>) :: l ->
                    if item.Class = lookingForClass then
                        let m' = (findAndIntegrate m item.Class item.GroupOfProteinIDs)
                        loop m' l remaining lookingForClass
                    else
                        loop m l (item::remaining) lookingForClass
            | [] ->
                if lookingForClass = PeptideEvidenceClass.C3b then
                    m
                else
                    loop m remaining [] (increaseClass lookingForClass)
        loop (ClassMap.empty()) proteinClassItems [] PeptideEvidenceClass.C1a

    let inferSequences integrationFunction (proteinClassItems: ProteinClassItem<'sequence> list) =
        ///Adds a value v to the set of values already associated with key k
        let internalAdd (dict:Dictionary<string,HashSet<'sequence>>) (k:string) (v:'sequence) =
            match dict.TryGetValue(k) with
            | (true, container) ->
                container.Add(v) |> ignore
            | (false,_) -> 
                let tmp = HashSet<'sequence>(new StructuralEqualityComparer<'sequence>())
                tmp.Add(v) |> ignore
                dict.Add(k,tmp) 
        let classes = groupIntoClasses integrationFunction proteinClassItems
        let seqDict = System.Collections.Generic.Dictionary<string,HashSet<'sequence>>()
        proteinClassItems
        |> List.iter (fun i -> internalAdd seqDict (proteinGroupToString i.GroupOfProteinIDs) i.PeptideSequence)
        classes
        |> Seq.collect (fun kv -> 
            let c = kv.Key
            kv.Value.GetArrayOfValues |> Seq.map (fun group -> createInferredProteinClassItem group c (Seq.toArray (seqDict.Item group)))
            )

            