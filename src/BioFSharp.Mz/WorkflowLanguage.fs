namespace BioFSharp.Mz

module WorkflowLanguage = 

    module Definition = 

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
    
        type IOType = 
            | Folder of System.IO.DirectoryInfo
            | File of System.IO.FileInfo
            | Files of System.IO.FileInfo list 
    
        type Operation<'a>  = {
            Id               : System.Guid
            Name             : string
            Operator         : string
            Input            : IOType
            Output           : IOType
            Parameters       : 'a       
            }
            
        let createProcessDescription id name operator input output parameters = { 
            Id = id; Name = name; Operator = operator; Input = input; Output = output; Parameters = parameters }
        
        type Workflow<'a> = {
            Processes : Operation<'a> list 
            }

        let createWorkflow processes = 
            {Processes = processes}

    module MSProcessing = 

        open Definition
        open Newtonsoft.Json

        type PeakPickingParams = {
            Ms1PeakPicking      : SignalDetection.Wavelet.WaveletParameters 
            PaddingParameters   : SignalDetection.Padding.PaddingParameters
            Ms2PeakPicking      : SignalDetection.Wavelet.WaveletParameters
            }
    
        let createPeakPickingParams ms1PeakPicking paddingParameters ms2PeakPicking = {
            Ms1PeakPicking      = ms1PeakPicking
            PaddingParameters   = paddingParameters
            Ms2PeakPicking      = ms2PeakPicking
            }        
        
        type MSParameters =
            | PeakPicking of PeakPickingParams

        ///
        let ToJSon (operation: Definition.Operation<MSParameters>) =
            Newtonsoft.Json.JsonConvert.SerializeObject operation
        
        ///
        let ofJson json = 
            Newtonsoft.Json.JsonConvert.DeserializeObject<MSParameters>(json)
