﻿namespace BioFSharp.Mz


type IPeak =    
    abstract member Mz        : float
    abstract member Intensity : float   


[<Struct>]
type Peak(mz:float,intensity:float) = 
        member this.Mz = mz
        member this.Intensity = intensity
        
        interface IPeak with
            member this.Mz = mz
            member this.Intensity = intensity
//
//
//type IonTypes =
//    | Unknown       =  0
//    | Precursor     =  1
//    | A             =  2
//    | B             =  3
//    | C             =  4
//    | X             =  5
//    | Y             =  6
//    | Z             =  7
//    | AlossH2O      =  8
//    | AlossNH3      =  9
//    | BlossH2O      = 10 
//    | BlossNH3      = 11 
//    | ClossH2O      = 12 
//    | ClossNH3      = 13 
//    | XlossH2O      = 14 
//    | XlossNH3      = 15 
//    | YlossH2O      = 16 
//    | YlossNH3      = 17
//    | ZlossH2O      = 18 
//    | ZlossNH3      = 19 
//    | Immonium      = 20 
//    | Neutral       = 21
//    | Diagnostic    = 22

type IonSeries =  
    | A 
    | B 
    | C 
    | X 
    | Y 
    | Z 

type IonTypes =
    | Unknown 
    | Precursor 
    | Immonium
    | Main         of IonSeries
    | LossH20      of IonSeries
    | LossNH3      of IonSeries
    | NeutralLoss  of IonSeries
    | Diagnostic   of IonSeries

type Tag<'t,'v> = {
    Meta : 't
    Data : 'v
    }


type PeakAnnotation = Tag<IonTypes,Peak>


module Peaks =
 
    ///
    let createPeak mzData intensityData = 
        let pk = new Peak(mzData, intensityData)
        pk
    
    ///
    let createTag meta data = {
        Meta=meta
        Data=data
        }
    
    ///
    let createPeakAnnotation (ionType:IonTypes) (peak:Peak) =
        let pA: PeakAnnotation = createTag ionType peak
        pA