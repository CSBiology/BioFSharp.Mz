namespace BioFSharp.Mz


module MSProcessing = 

    open BioFSharp.WorkflowLanguage 
    open Newtonsoft.Json

    type DataFormat = 
        | Wiff 
        | Baf 
        | ThermoRawFile
        | CSV
        | XML 
        | SQLite 
        | Txt
        | Custom 

    type DataModel =  
        | Wiff of DataFormat 
        | Baf of DataFormat 
        | ThermoRawFile of DataFormat 
        | MzML of DataFormat 

    type PeakPickingParams = {
        CompressData        : bool
        Ms1PeakPicking      : SignalDetection.Wavelet.WaveletParameters option
        PaddingParameters   : SignalDetection.Padding.PaddingParameters option
        Ms2PeakPicking      : SignalDetection.Wavelet.WaveletParameters option 
        }
    
    let createPeakPickingParams compressData ms1PeakPicking paddingParameters ms2PeakPicking = { 
        CompressData = compressData
        Ms1PeakPicking      = ms1PeakPicking
        PaddingParameters   = paddingParameters
        Ms2PeakPicking      = ms2PeakPicking
        }        
        
    type MSParameters =
        | PeakPicking of PeakPickingParams
            
    ///
    let operationToJSon (operation: Definition.Operation<MSParameters>) =
        Newtonsoft.Json.JsonConvert.SerializeObject operation
        
    ///
    let operationOfJson json = 
        Newtonsoft.Json.JsonConvert.DeserializeObject<Definition.Operation<MSParameters>>(json)

    ///
    let workFlowToJSon (operation: seq<Definition.Operation<MSParameters>>) =
        Newtonsoft.Json.JsonConvert.SerializeObject operation
        
    ///
    let workFlowToJson json = 
        Newtonsoft.Json.JsonConvert.DeserializeObject<seq<Definition.Operation<MSParameters>>>(json)

        