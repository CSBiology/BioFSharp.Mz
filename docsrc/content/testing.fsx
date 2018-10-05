(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
//#I "../../bin"
//TODO: Add FILES and LINKS to other tutorials
//#I "../../bin/BioFSharp.Mz/net47/"
#I @"C:\Users\david\Source\Repos\netCoreRepos\BioFSharp.Mz\src\BioFSharp.Mz\bin\Release\net47/"
#r "BioFSharp.dll"
#r "BioFSharp.Mz.dll"
#r "FSharpAux.dll"
#r "FSharp.Stats.dll"
#r "Newtonsoft.Json.dll"
#r @"C:\Users\david\Source\Repos\netCoreRepos\BioFSharp.Mz\packages\FSharp.Plotly\lib\net47\FSharp.Plotly.dll"
(**
Peptide retrieval by mass
=========================
*)

(**
This part of the documentation aims to give you a brief overview of the functionality of SearchDB.fs. This Module contains functions to initialize a user configured peptide database and to subsequently perform common database operations such as inserts or lookUps.


Creating a peptide database
---------------------------

Prior to database creation the user has to specify the database parameters. The following examples will show this process in detail. 
Some complex parameters will be discussed more specifically.

The concept of peptide modifications
------------------------------------
In the library the term "Modification" is understood as literally any change to an amino acid formula. These changes can vary 
from the exchange of an isotope (e.g. a N14 gets replaced with a N15) to an addition/deletion of chemical groups e.g. a phosphorylation.

Differentiation of fixed, variable and isotopic modifications
-----------------------------------------------------------
Creating a peptide database the user not only has to specify the modification itself but also its persistance. This is due to the fact 
that the abundance of a modification can vary from one sample to another. This thought lead to a distinction between fixed, variable and 
isotopic modifications during the construction of this library. 

In a natural sample it is very unlikely for a modification like a phosphorylation to be present at every matching amino acid of every
peptide. The conclusion of this finding is that a database has to contain both, a modified and a unmodified version of the peptide simulating a 
**variable** occurency of the modification coining the term **"VariableModification"**.   

Contrary, further manipulations to the sample e.g. a selection of only N-terminally phosphorylated peptides can result in a sample
for which it can be assumed that roughly every peptide contains this modification. This would justify a creation of a peptide database with a 
**fixed** N-terminal phosphorylation at every matching amino acid of every peptide. Hereby coining the term **"FixedModification"**.

-----------------------------------------------------------
This expression constructs a new instance of the type SearchModification, which is needed to define Fixed- and VariableModifications
This type is inferred of the XMod ontology and the information stored in this type is needed to build instances of the type Modification.     
*)
/// computes a non-circular convolution of x and y 
//let convolve (x:'a []) (y:'a []):'a []= 
//    let k = (x.Length+y.Length-1)
//    let tmp = Array.zeroCreate k
//    for i = 0 to x.Length-1 do
//        for j = 0 to (y.Length-1) do 
//            tmp.[i+j] <- tmp.[i+j] + (x.[i] * y.[j])
//    tmp 

//let x = [|0. ..100.|]
"C:\Users\david\OneDrive - tukl\Dokumente\BioInfo\82_TorstenProteins\CsbScaffold\.env\packages\runtime.opensuse.13.2-x64.runtime.native.System.Security.Cryptography.OpenSsl\runtime.opensuse.13.2-x64.runtime.native.system.security.cryptography.openssl.4.3.3.nupkg".Length
//let y = x |> Array.map (FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue (FSharp.Stats.Vector.ofArray [|100.;50.;4.|]))
//let discFilter = [|0.;0.;0.;0.;1.;-1.;0.;0.;0.;0.;0.|]
//let yy = 
//    let tmp =
//        (convolve y discFilter)
//    tmp.[5..tmp.Length-6]
//x.Length
//yy.Length

//open FSharp.Plotly
//[
//Chart.Point(x,yy)
//Chart.Point(x,y)
//]
//|> Chart.Combine
//|> Chart.Show



open BioFSharp
open BioFSharp.Mz
open ModificationInfo
open AminoAcids
open SearchDB
open BioFSharp.Digestion
#time

open System 
open System.Runtime.Serialization.Formatters.Binary
open System.IO
open BioFSharp.Mz.SearchDB




/// Returns a instance of the type SearchModification
let phosphorylation = {
    // Name of this SearchModification 
    Name="Phosphorylation" 
    //
    IsBiological = true
    // Xmod accession Number
    Accession="21" 
    // A brief description of the SearchModification
    Description="Addition of a phosphate group to the AA residue"
    // A string representation of the formula composition 
    Composition= "HO3P"
    // Specifices the ModLocation
    Site= [Specific(Met,ModLocation.Residual);(* Specific(Thr,ModLocation.Residual); Specific(Tyr,ModLocation.Residual)*)]
    // Specifies if the Modification is additive or subtractive
    MType=SearchModType.Plus
    // If a modified AminoAcid is shown via string representation this code is placed before
    // the AminoAcid Symbol e.g. "phS"
    XModCode= "ph"
    } 

(**

As mentioned before, there is a need to model modifications that affect the peptide sequences on a **isotopic** level. This could be understood 
as a special kind of FixedModification with the additional constrain that there is no addition or deletion of a chemical group but a isotopic 
modification affecting every amino acid. e.g. a exchange of every N14 to a N15. As one could expect, these features coined the term "IsotopicModification".   

This expression returns a instance of the type SearchInfoIsotopic, which is needed to define IsotopicModifications. This type is needed if the user 
plans to introduce a IsotopicModifications e.g a replacement of the N14 isotopes with heavier N15 isotopes.
*)

/// Returns a instance of the type SearchInfoIsotopic 
let isoMod = (SearchDB.createSearchInfoIsotopic "N15" Elements.Table.N Elements.Table.Heavy.N15)

let isoM = isoMod |> SearchDB.createIsotopicMod
let sample = BioList.ofAminoAcidString "AAA"
let isoModSample = sample |> List.map (fun x -> AminoAcids.setModification isoM x)

sample.[0] = isoModSample.[1]

let getFrags seque = 
    (Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList Formula.monoisoMass seque).TargetMasses 
    |> List.map (fun x -> x.MainPeak)
getFrags sample

for i = 1 to 10000 do
    getFrags isoModSample

(**
Selecting and initializing a mass function
------------------------------------------
The following expression returns a instance of a massfunction. The user can choose between functions that either return 
the monoisotopic or the average mass of types that inherit the IBioItem interface (e.g. AminoAcids).  
Moreover, each function exists in a "WhitMemP" version that uses a dictionary to "remember" (store) already computed masses. 
This approach is especially useful when the user aims to do many mass computations subsequently, as it is the case when
building a new database. The usage of a dictionary is the reason for this function to be initialized prior to database
creation.      
*)

///Returns a instance of the massfunction BioItem.initMonoisoMassWithMemP
let massF : (IBioItem -> float) =  BioItem.initMonoisoMassWithMemP


//let sample = BioList.ofAminoAcidString "AAA"

let m = AminoAcid.Met

let acetylation'NTerm' =
    createSearchModification  "Acetyl(N-Term)"  "1" "Acetylation of the protein N-terminus" false "C2H2ON"
        [Any(ModLocation.Nterm)] SearchModType.Plus "no"

let acetylation'NTermbio' =
    createSearchModification  "Acetyl(N-Termbio)"  "2" "Acetylation of the protein N-terminusss" true "C2H2ON"
        [Any(ModLocation.Nterm)] SearchModType.Plus "ja"

let nonBioOx = SearchDB.ModCombinator.convertSearchModification [AminoAcid.Met] acetylation'NTerm' |> List.item 0 |> snd
let bioOx    = SearchDB.ModCombinator.convertSearchModification [AminoAcid.Met] acetylation'NTermbio' |> List.item 0 |> snd


let mIso = AminoAcid.Met |> AminoAcids.setModification isoM

let mIsoBioOx    = mIso |> AminoAcids.setModification bioOx
let mIsononBioOx = mIso |> AminoAcids.setModification nonBioOx

//for i = 1 to 1000000 do 
//let modW = "C2H2ON" |> Formula.parseFormulaString |> Formula.monoisoMass
//massF m 
//|> (+) modW 

//massF mIso
//massF mIsononBioOx
//massF mIsoBioOx

//massF isoModSample.[1]

///
let calcAIonFragMass (massfunction:Formula.Formula -> float) acc aa =
    acc + massfunction (AminoAcids.formula aa) - (massfunction (AminoAcids.isotopicLabelFunc aa Formula.Table.CO)) 
    
let calcAIonFragMassLol (massfunction:IBioItem -> float) acc aa =
    acc + (massfunction (aa |> AminoAcids.setModification bioOx))
         
calcAIonFragMass Formula.monoisoMass 0. mIsononBioOx

//for i = 1 to 1000000 do 
//    calcAIonFragMassLol massF 0. mIsononBioOx
//calcAIonFragMass Formula.monoisoMass 0. mIsononBioOx

(**
Incorporation of all parameters in a single record type
-------------------------------------------------------
This type contains all formerly designed Parameters and additional parameters. 
*)

let mass:IBioItem -> float =  BioItem.monoisoMass //BioItem.initMonoisoMassWithMemP

let moddddd = GlobalModificationInfo.initGlobalModificationDeltaOfMod mass [SearchDB.createIsotopicMod isoMod]

let acetylation'ProtNTerm' =
    createSearchModification "Acetyl(Protein N-Term)" "1" "Acetylation of the protein N-terminus" true "C2H2O"
        [Any(ModLocation.ProteinNterm)] SearchModType.Plus "pn"


let acetylation'CTerm' =
    createSearchModification "Acetyl(C-Term)" "1" "Acetylation of the protein N-terminus" true "C2H2O"
        [Any(ModLocation.Cterm)] SearchModType.Plus "ac"

let acetylation'ProtCTerm' =
    createSearchModification "Acetyl(Protein C-Term)" "1" "Acetylation of the protein N-terminus" true "C2H2O"
        [Any(ModLocation.ProteinCterm)] SearchModType.Plus "pc"

/// Returns a instance of the type SearchDbParams
let paramTestN15 = {
        Name="Creinhardtii_fullN15_ProteinNTerm_final2_biofinal"
        // Path of db storage folder
        DbFolder            = (__SOURCE_DIRECTORY__ + "/data/")
        // Path of the user selected fasta file  
        FastaPath           = (__SOURCE_DIRECTORY__ + "/data/Chlamy_JGI5_5(Cp_Mp).fasta")
        // Function that specifies the conversion of the Fasta Header to the protein name inserted in the database
        FastaHeaderToName   = id
        // Protease that is used during in silico digestion of the proteins
        Protease            = Digestion.Table.Trypsin
        // Minimum number of missed cleavages
        MinMissedCleavages  = 0
        // Maximum number of missed cleavages
        MaxMissedCleavages  = 2
        // Maximum peptide mass inserted in the database
        MaxMass             = 15000.
        // Mininmum peptide length 
        MinPepLength        = 4
        // Maximum pepide length
        MaxPepLength        = 65
        // Contains a instance of SearchInfoIsotopic, set to None if no isotopic modification shall be introduced
        IsotopicMod           = [isoMod] 
        // MassMode used during db creation
        MassMode            = SearchDB.MassMode.Monoisotopic
        // Massfunction used during db creation
        MassFunction        = massF
        // List of "Fixedmodifications". 
        FixedMods           = [(*acetylation'NTerm';acetylation'NTermbio'*)]
        // List of "VariableModificaions".  
        VariableMods        = [Table.oxidation'Met';Table.acetylation'ProtNTerm'(*BioFSharp.Mz.SearchDB.Table.oxidation'Met';(*phosphorylation;*)acetylation'ProtNTerm'*)]
        // Maximum number of variable modifications per peptide this constrain is used to limit the computational 
        VarModThreshold     = 2
        }
#time  

(**

The following example illustrates how a database lookUp can be initialized. The function SearchDB.getPeptideLookUpBy 
takes the defined SearchDbParams and returns a function that takes two float values (representing masses) as input 
parameters and returns a list of peptides each with a mass in the range of first to the second mass used as input 
parameters. 

SearchDB.getPeptideLookUpBy first checks if a database with the given parameters already exists before it starts the 
database creation process. If there is a database matching the user given SearchParameters the process will be finished 
nearly immediatly. If a new database has to be created the computation time will depend on the complexity of the given 
parameters such as the size of the FASTA file, the number of given Modifications, the maxPepLength et cetera. 
5+5   
*)
/// Returns a function that takes two masses as input parameters and returns a list of peptides (wrapped in the type LookUpResult)
let n15DBLookUpBy =  
        SearchDB.getPeptideLookUpBy paramTestN15   

n15DBLookUpBy 1000. 1001.

let Params = 
    SearchDB.Db.isExistsBy paramTestN15


#r "System.Data.SQLite.dll"
open System.Data.SQLite
/// Returns true if a db exists with the same parameter content
let isExistsBy (sdbParams:SearchDbParams) =       
    let fileName = SearchDB.Db.getNameOf sdbParams
    let connectionString = sprintf "Data Source=%s;Version=3" fileName
    use cn = new SQLiteConnection(connectionString)
    cn.Open()
    SearchDB.Db.selectSdbParamsby cn sdbParams 

let Params' = 
    paramTestN15

Params'.FixedMods@Params'.VariableMods      


BioItemsConverter.OptionConverter.charToOptionAminoAcid 'c'
open BioItemsConverter
open BioList
open SearchDB
open System
open FSharp.Stats.SpecialFunctions

Char.IsUpper 'a'
     

let stringF = SearchDB.initOfModAminoAcidString paramTestN15.IsotopicMod (paramTestN15.FixedMods@paramTestN15.VariableMods)

let res = stringF 0 "RISLKPKF*"
let resIso = stringF 1 "RISLKPKF*"
res |> BioSeq.toString
BioSeq.toMonoisotopicMassWith (Formula.Table.H2O |> Formula.monoisoMass) res
BioSeq.toMonoisotopicMassWith (Formula.Table.H2O |> Formula.monoisoMass) resIso

//    /// Generates amino acid sequence of one-letter-code string containing modified AminoAcids indicated by 2 lowercase digits per modification. 
//let ofModAminoAcidString (converter: OptionConverter.AminoAcidOptionConverter) (converter': 'b -> Modification ) (isotopMod: Modification list option) (xModToSearchMod: Map<string,'b>) (aaStr: string) : BioList<_>  =
//    let aaStrL = aaStr.Length
//    let converter a = 
//        match isotopMod with 
//        | Some modL -> converter a |> Option.get |> (setModifications modL)
//        | None      -> converter a |> Option.get
//    let accumMods startIdx =
//        let rec loop currentIdx acc (sequence: string) =
//            if currentIdx < 0 then acc
//            elif sequence.[currentIdx] = '[' then 
//                match Map.tryFind acc xModToSearchMod with 
//                | Some x -> x
//                | None   ->  
//            else 
//                loop (currentIdx-1) (string sequence.[currentIdx] + acc ) sequence 
//        loop startIdx "" sequence
//    let rec loopWithGlobal count (modAcc: 'b list) acc (converter:OptionConverter.AminoAcidOptionConverter) (xModToSearchMod: Map<string,'b>)  (aaStr: string) = 
//        if count = aaStrL then 
//                acc
//        else 
//                let currentC = aaStr.[count]
//                if  (currentC |> Char.IsUpper = true) && modAcc = [] then 
//                    let currentA = (converter currentC).Value 
//                                |> setModifications isotopMod.Value
//                    loopWithGlobal (count+1) [] (currentA::acc) converter xModToSearchMod aaStr 
//                elif 
//                    ((currentC |> Char.IsUpper) = true) then
//                    let modList = List.map converter' (modAcc)
//                    let tmpAa = setModifications (isotopMod.Value@modList) (converter currentC).Value
//                    loopWithGlobal (count+1) [] (tmpAa::acc) converter xModToSearchMod aaStr
//                else 
//                    match Map.tryFind aaStr.[count.. count+1] xModToSearchMod with
//                    | Some modi -> loopWithGlobal (count+1) (modi::modAcc) acc converter xModToSearchMod aaStr
//                    | None      -> loopWithGlobal (count+1) modAcc acc converter xModToSearchMod aaStr                 
//    match isotopMod.IsSome with
//    | true  -> loopWithGlobal 0 [] [] converter xModToSearchMod aaStr
//    | false -> ofRevModAminoAcidString converter converter' xModToSearchMod aaStr



    ///// Generates amino acid sequence of one-letter-code string containing modified AminoAcids indicated by 2 lowercase digits per modification. 
    //let ofRevModAminoAcidString (converter: OptionConverter.AminoAcidOptionConverter) (converter': 'b -> Modification ) (xModToSearchMod: Map<string,'b>) (aaStr: string) : BioList<_>  =
    //    let aaStrL = aaStr.Length
    //    let rec loop count (modAcc: 'b list) acc (converter:OptionConverter.AminoAcidOptionConverter) (xModToSearchMod: Map<string,'b>)  (aaStr: string) = 
    //        if count = aaStrL then 
    //             acc
    //        else 
    //             let currentC = aaStr.[count]
    //             if  (currentC |> Char.IsUpper = true) && modAcc = [] then 
    //                 loop (count+1) modAcc ((converter currentC).Value::acc) converter xModToSearchMod aaStr 
    //             elif 
    //                 ((currentC |> Char.IsUpper) = true) then
    //                 let modList =
    //                       List.map converter' modAcc
    //                 let tmpAa = setModifications modList (converter currentC).Value
    //                 loop (count+1) [] (tmpAa::acc) converter xModToSearchMod aaStr
    //             else 
    //                 match Map.tryFind aaStr.[count.. count+1] xModToSearchMod with
    //                 | Some modi -> loop (count+1) (modi::modAcc) acc converter xModToSearchMod aaStr
    //                 | None      -> loop (count+1) modAcc acc converter xModToSearchMod aaStr 
        
    //    loop 0 [] [] converter xModToSearchMod aaStr


    ///// Generates amino acid sequence of one-letter-code string containing modified AminoAcids indicated by 2 lowercase digits per modification. 
    //let ofRevModAminoAcidStringWithIsoMod (converter: OptionConverter.AminoAcidOptionConverter) (converter': 'b -> Modification ) (isotopMod: Modification list option) (xModToSearchMod: Map<string,'b>) (aaStr: string) : BioList<_>  =
    //    let aaStrL = aaStr.Length
    //    let rec loopWithGlobal count (modAcc: 'b list) acc (converter:OptionConverter.AminoAcidOptionConverter) (xModToSearchMod: Map<string,'b>)  (aaStr: string) = 
    //        if count = aaStrL then 
    //             acc
    //        else 
    //             let currentC = aaStr.[count]
    //             if  (currentC |> Char.IsUpper = true) && modAcc = [] then 
    //                 let currentA = (converter currentC).Value 
    //                                |> setModifications isotopMod.Value
    //                 loopWithGlobal (count+1) [] (currentA::acc) converter xModToSearchMod aaStr 
    //             elif 
    //                 ((currentC |> Char.IsUpper) = true) then
    //                 let modList = List.map converter' (modAcc)
    //                 let tmpAa = setModifications (isotopMod.Value@modList) (converter currentC).Value
    //                 loopWithGlobal (count+1) [] (tmpAa::acc) converter xModToSearchMod aaStr
    //             else 
    //                 match Map.tryFind aaStr.[count.. count+1] xModToSearchMod with
    //                 | Some modi -> loopWithGlobal (count+1) (modi::modAcc) acc converter xModToSearchMod aaStr
    //                 | None      -> loopWithGlobal (count+1) modAcc acc converter xModToSearchMod aaStr                 
    //    match isotopMod.IsSome with
    //    | true  -> loopWithGlobal 0 [] [] converter xModToSearchMod aaStr
    //    | false -> ofRevModAminoAcidString converter converter' xModToSearchMod aaStr

(**
This example shows how the pattern used before can be modified slightly to return a function that only needs one parameter as 
input to return a LookUpResult list. 
*)

/// Returns a function that takes one inputMass and returns a list of peptides (wrapped in the type LookUpResult) in the range 
/// of the inputMass-0.5 to inputMass+0.5 
let n15LookUpPeptideBy inputMass =
        let lowerBorder = inputMass-0.5
        let upperBorder = inputMass+0.5 
        n15DBLookUpBy lowerBorder upperBorder

n15LookUpPeptideBy 1000.1
(**
This function can be used as shown in the code snippet below. As a inputMass the monoisotopic mass of the peptide *ANLGMEVMHER* can be used as a input. 
This peptide is encoded by a protein which can be found in the FASTA used for the database creation. Surprisingly the list you obtain will be empty because 
the database does not contain any peptides in the given massRange.   
*)

/// PeptideMass of the peptide known from the peptide *ANLGMEVMHER*. 
let peptideMass = 1285.590726

(**
The possibility that the database is empty can be ruled out by performing a lookUp with the mass 1000., which should return a list filled with 
items of the type LookUpResult. The explanation to this behaviour can be found examining the SearchDbParameters bound to the name **paramTestN15**.
When looking at the parameters we see that we introduced a IsotopicModification. If the peptide was optained from a sample without this isotopic 
modification, a second lookUp in a database lacking IsotopicModification might return a result list containing the peptide of interest.

As a first step a new instance of the type SearchDbParams has to be instantiated and a database has to be created:
*)

/// Returns a instance of the type SearchDbParams
let paramTestN14 = {
        Name="Creinhardtii_N1422"
        DbFolder            = (__SOURCE_DIRECTORY__ + "/data/")
        FastaPath           = (__SOURCE_DIRECTORY__ + "/data/Chlamy_JGI5_5(Cp_Mp)_truncated.fasta")
        FastaHeaderToName   = id
        Protease            = Digestion.Table.getProteaseBy "Trypsin"
        MinMissedCleavages  = 0
        MaxMissedCleavages  = 3
        MaxMass             = 15000.
        MinPepLength        = 4
        MaxPepLength        = 65
        IsotopicMod           = []
        MassMode            = SearchDB.MassMode.Monoisotopic
        MassFunction        = massF
        FixedMods           = [] 
        VariableMods        = [BioFSharp.Mz.SearchDB.Table.oxidation'Met']
        VarModThreshold     = 1
        }
#time
/// Returns a function that takes two masses as input parameters and returns a list of peptides (wrapped in the type LookUpResult)
let n14DBLookUpBy =  
        SearchDB.getPeptideLookUpBy paramTestN14   
/// Returns a function that takes one inputMass and returns a list of peptides (wrapped in the type LookUpResult) in the range 
/// of the inputMass-0.5 to inputMass+0.5 
let n14LookUpPeptideBy inputMass =
        let lowerBorder = inputMass-0.5
        let upperBorder = inputMass+0.5 
        n14DBLookUpBy lowerBorder upperBorder

(**
After these steps the following lookup will result in list of lookupResults containing the peptide "ANLGMEVMHER".
*)

// Returns the result of the peptide lookup
let n14LookUpResult = n14LookUpPeptideBy peptideMass