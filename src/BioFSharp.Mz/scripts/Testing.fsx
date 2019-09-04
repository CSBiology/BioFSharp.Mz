#r "netstandard"
#r @"../../../\packages\System.Data.SQLite.Core\lib\net46\System.Data.SQLite.dll"
#r @"../../../packages\BioFSharp\lib\netstandard2.0\BioFSharp.dll"
#r @"../../../bin\BioFSharp.Mz\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.dll"
open BioFSharp
open BioFSharp.Mz

let massF:(IBioItem->float) = BioItem.initMonoisoMassWithMem

let BioSeq = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString
let BioSeqH = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString |> List.map (fun x ->AminoAcids.setModifications [ModificationInfo.Table.N15] x)

BioSeqH |> BioList.toFormula |> Formula.toString
let mono = BioSeqH |> BioList.toMonoisotopicMassWith (Formula.Table.H2O |> Formula.monoisoMass)
ModificationInfo.Table.NH3_loss |> massF

let res  = Fragmentation.Series.bOfBioList massF BioSeq
let res2 = Fragmentation.Series.bOfBioList massF BioSeqH
let light = res.[0].MainPeak.Mass
let lightLoss = res.[0].DependentPeaks.[0].Mass
light - lightLoss 

let heavy = res2.[0].MainPeak.Mass

let deph20Loss = res2.[0].DependentPeaks.[0].Mass
open Elements
System.Numerics.Complex(-0.0926829268292669, 22.0592593273255)

open System

let compareComplexNumbers (c1:System.Numerics.Complex) (c2:System.Numerics.Complex) = 
    (Math.Abs(c1.Real - c2.Real) <= System.Double.Epsilon) &&
    (Math.Abs(c1.Imaginary - c2.Imaginary) <= System.Double.Epsilon);
      
let rootO = 
    match Elements.Table.O with 
    | Tri o -> o.Root
    //| _ -> 




let rootO' =
    System.Numerics.Complex(-0.0926829268292669, 22.0592593273255),
    System.Numerics.Complex(-0.0926829268292696, -22.0592593273255)
let Zn = Multi (createMulti "Zn" (Isotopes.Table.Zn64,Isotopes.Table.Zn64.NatAbundance ) (Isotopes.Table.Zn66,Isotopes.Table.Zn66.NatAbundance) (Isotopes.Table.Zn68,Isotopes.Table.Zn68.NatAbundance) [|Isotopes.Table.Zn67;Isotopes.Table.Zn70|])
open Formula
"H2O" |> parseFormulaString |> Formula.monoisoMass
open Formula
compareComplexNumbers (fst rootO) (fst rootO)
//let iso = getSinglePhiM Table.Na 1000000. 2.
//let iso = getSinglePhiM Table.H 1000000. 2.
//let iso = getSinglePhiM Table.O 1000000. 2.
//let iso = getSinglePhiM Zn 1000. 3.
open FSharpAux
let massF : IBioItem -> float= BioItem.initMonoisoMassWithMemP  
[0.. 1000]
|> List.map (fun x -> massF AminoAcids.Ala)
[0.. 10000000]
|> PSeq.map (fun x -> massF AminoAcids.Ala)
|> Array.ofSeq
open IsotopicDistribution.MIDA
/// normalizes istopic distribution probabilities to sum up to 1.
let normalizeByProbSum minProb (pols: Polynomial list) = 
    let sum = 
        pols
        |> List.filter (fun x -> x.Probability > minProb)
        |> List.sumBy (fun x -> x.Probability)
    pols 
    |> List.map (fun x -> {x with Probability = x.Probability/sum})


"LTYYTPDYVVR"
|> BioList.ofAminoAcidString
|> BioList.toFormula
//|> fun f -> Formula.lableElement f Elements.Table.N Elements.Table.Heavy.N15
//|> (IsotopicDistribution.MIDA.ofFormula (normalizeByProbSum) 0.01 (exp(-20.)) 2)
|> (IsotopicDistribution.BRAIN.ofFormula 10 )



"LTYYTPDYVVR"
|> BioList.ofAminoAcidString
|> BioList.toFormula
|> (IsotopicDistribution.BRAIN.ofFormula 25)

BioItem.monoisoMass AminoAcids.Ala 
   
   
//[ Array.iter (fun x -> if x > 0. then do yield x ) [|-50. .. 50.|]]


let f = "N100C12O8" |> parseFormulaString 
#time


open System.Data.SQLite
let cn = SearchDB.getDBConnection @"C:\Users\david\source\repos\ProteomIQon\src\ProteomIQon\Scripts\GiadaTest.db"
let p = SearchDB.getSDBParamsBy @"C:\Users\david\source\repos\ProteomIQon\src\ProteomIQon\Scripts\GiadaTest.db"
let look = SearchDB.getThreadSafePeptideLookUpFromFileBy cn p

let psm = look 1001. 1002.

let modPep = 
    psm
    |> List.find (fun x -> x.GlobalMod = 0)

modPep.Mass
modPep.BioSequence |> List.sumBy massF |> (+) (Formula.monoisoMass Formula.Table.H2O)  

let resH = Fragmentation.Series.bOfBioList massF modPep.BioSequence 
let light' = resH.[0].MainPeak.Mass
let Loss' = resH.[0].DependentPeaks.[0].Mass

light' - Loss' 


let heavy = res2.[0].MainPeak.Mass

let deph20Loss = res2.[0].DependentPeaks.[0].Mass
heavy - deph20Loss  



//let light = res.[0].MainPeak.Mass
let heavywaterl = res2.[0].DependentPeaks.[0].Mass

heavy - light
heavywaterl - heavy


let ala = AminoAcids.Ala |> massF
let aminoLoss = 
    AminoAcids.setModification ModificationInfo.Table.NH3_loss AminoAcids.Ala
    |> massF

ala - aminoLoss 

//let alaH = 
//for i = 10 to 100000 do 
//        //AminoAcids.setModification ModificationInfo.Table.N15 AminoAcids.Ala
//        AminoAcids.Ala |> massF 
    //ModificationInfo.Table.NH3_loss 
    //|> massF

   
////let aminoLossH = 
//{ModificationInfo.Table.NH3_loss with Name= "Haralddddd"; Modify =  ModificationInfo.Table.NH3_loss.Modify >> ModificationInfo.Table.N15.Modify } 
//|> massF

//let x = {ModificationInfo.Table.NH3_loss with Name="moped"; Modify =  >>  ModificationInfo.Table.N15 .Modify } 
//x|> massF


//AminoAcids.Ala = (AminoAcids.Ala |> (AminoAcids.setModification ModificationInfo.Table.NH3_loss)  )
//(AminoAcids.Ala |> (AminoAcids.setModification ModificationInfo.Table.NH3_loss)  ) = (AminoAcids.Ala |> (AminoAcids.setModification ModificationInfo.Table.NH3_loss)  )
//(AminoAcids.Ala |> (AminoAcids.setModification ModificationInfo.Table.NH3_loss)  ) = (AminoAcids.Ala |> (AminoAcids.setModification {ModificationInfo.Table.NH3_loss with Name="moped"; Modify = ModificationInfo.Table.NH3_loss.Modify >>  ModificationInfo.Table.N15 .Modify })  )


//// Modifies existing modification (md) using (md')
//let modifyModification (md':ModificationInfo.Modification) (md:ModificationInfo.Modification) =
//    {md with Name = sprintf "%s_%s" md.Name md'.Name ; Modify = md.Modify >>  md'.Modify} 



//let heavyNh3 = modifyModification ModificationInfo.Table.N15 ModificationInfo.Table.NH3_loss

//ModificationInfo.Table.NH3_loss |> massF

//heavyNh3 |> massF

//ModificationInfo.Table.N15 |> massF




