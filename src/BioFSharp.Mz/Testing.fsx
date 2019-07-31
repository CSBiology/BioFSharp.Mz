#r "netstandard"

#r @"../../packages\BioFSharp\lib\netstandard2.0\BioFSharp.dll"
#r @"../../bin\BioFSharp.Mz\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.dll"
open BioFSharp
open BioFSharp.Mz

let massF:(IBioItem->float) = BioItem.initMonoisoMassWithMem

let BioSeq = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString
let BioSeqH = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString |> List.map (fun x ->AminoAcids.setModifications [ModificationInfo.Table.N15] x)
BioSeqH |> BioList.toFormula |> Formula.toString
let mono = BioSeqH |> BioList.toMonoisotopicMassWith (Formula.Table.H2O |> Formula.monoisoMass)
ModificationInfo.Table.NH3_loss |> massF

let res  = Fragmentation.Series.yOfBioList massF BioSeq
let res2 = Fragmentation.Series.yOfBioList massF BioSeqH

let light = res.[0].MainPeak.Mass
let heavy = res2.[0].MainPeak.Mass

let light = res.[0].MainPeak.Mass
let heavywaterl = res2.[0].DependentPeaks.[1].Mass

heavy - light
heavywaterl - heavy