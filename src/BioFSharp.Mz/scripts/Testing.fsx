#r "netstandard"

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




