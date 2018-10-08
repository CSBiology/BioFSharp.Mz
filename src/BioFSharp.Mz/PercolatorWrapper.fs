﻿namespace BioFSharp.Mz

module PercolatorWrapper = 

    module Parameters = 

        // TODO: Refactor in FSharp.Care 
        let fileInfoToLinuxPath (fI:System.IO.FileInfo) =
            let directory = fI.Directory.Root.ToString().Substring(0,1).ToLower()
            "/mnt/" + directory + (fI.FullName.Substring(2).Replace('\\', '/'))
        
        let fileInfoToWindowsPath (fI:System.IO.FileInfo) =
            fI.FullName

        type GeneralOptions = 
            /// Display the help message
            | Help
            /// Set verbosity of output: 0=no processing info, 5=all. Default = 2
            | VerbosityOfOutput of int
            /// Do not remove redundant peptides, keep all PSMs and exclude peptide level probabilities.
            | OnlyPSMs
            /// Use the mix-max method to assign q-values and PEPs. Note that this option only has an 
            /// effect if the input PSMs are from separate target and decoy searches. This is the default setting.
            | PostProcessing_MIXMAX
            ///Replace the mix-max method by target-decoy competition for assigning q-values and PEPs. 
            /// If the input PSMs are from separate target and decoy searches, Percolator's SVM scores 
            /// will be used to eliminate the lower scoring target or decoy PSM(s) of each scan+expMass combination. 
            /// If the input PSMs are detected to be coming from a concatenated search, this option will be turned on automatically, 
            /// as this is incompatible with the mix-max method. In case this detection fails, turn this option on explicitly.
            | PostProcessing_TargetDecoyCompetition

        let private stringOfgO (gO:GeneralOptions) =
            match gO with 
            | Help                                   -> " --help" 
            | VerbosityOfOutput nr                   -> " --verbose " + string nr   
            | OnlyPSMs                               -> " --only-psms " 
            | PostProcessing_MIXMAX                  -> " --post-processing-mix-max "
            | PostProcessing_TargetDecoyCompetition  -> " --post-processing-tdc "

        type FileInputOptions = 
            /// Read percolator tab-input format (pin-tab) from standard input.
            | PINTAB  of System.IO.FileInfo 
            /// Read percolator xml-input format (pin-xml) from standard input.
            | PINXML  of System.IO.FileInfo 
            /// Input file given in deprecated pin-xml format generated by e.g. sqt2pin with the -k option.
            | DeprecatedPINXML of System.IO.FileInfo
            /// Skip validation of pin-xml input file against xml schema.
            | SkipSchemeValidation
        
        let private stringOfFIO fileInfoToPath (fIO:FileInputOptions) =
            match fIO with 
            | PINTAB                info -> fileInfoToPath info 
            | PINXML                info -> fileInfoToPath info 
            | DeprecatedPINXML      info -> " --xml-in " + fileInfoToPath info 
            | SkipSchemeValidation ->       " --no-schema-validation "
    
        type FileOutputOptions =
            /// Output tab delimited results of peptides to a file instead of stdout (will be ignored if used with -U option)
            | POUTTAB_Peptides of System.IO.FileInfo
            /// Output tab delimited results for decoy peptides into a file (will be ignored if used with -U option).
            | POUTTAB_DecoyPeptides of System.IO.FileInfo
            /// Output tab delimited results of PSMs to a file instead of stdout.
            | POUTTAB_PSMs of System.IO.FileInfo
            /// Output tab delimited results for decoy PSMs into a file
            | POUTTAB_DecoyPSMs of System.IO.FileInfo
            /// Output tab delimited results of proteins to a file instead of stdout (Only valid if option -A or -f is active)
            | POUTTAB_Proteins of System.IO.FileInfo
            /// Output tab delimited results for decoy proteins into a file (Only valid if option -A or -f is active)
            | POUTTAB_DecoyProteins of System.IO.FileInfo
            /// Output computed features to given file in pin-tab format. Can be used to convert pin-xml to pin-tab.
            | POUTTAB_Features of System.IO.FileInfo
            /// Path to xml-output (pout) file.
            | POUTXML of System.IO.FileInfo
            /// Include decoys (PSMs, peptides and/or proteins) in the xml-output. Only available if -X is set.
            | IncludeDecoysInXML 

        let private stringOfFOO fileInfoToPath (fOO: FileOutputOptions) = 
            match fOO with 
            | POUTTAB_Peptides info      -> " --results-peptides " + fileInfoToPath info 
            | POUTTAB_DecoyPeptides info -> " --decoy-results-peptides " + fileInfoToPath info 
            | POUTTAB_PSMs info          -> " --results-psms " + fileInfoToPath info 
            | POUTTAB_DecoyPSMs info     -> " --decoy-results-psms " + fileInfoToPath info 
            | POUTTAB_Proteins info      -> " --results-proteins " + fileInfoToPath info 
            | POUTTAB_DecoyProteins info -> " --decoy-results-proteins " + fileInfoToPath info 
            | POUTTAB_Features info      -> " --tab-out " + fileInfoToPath info 
            | POUTXML info               -> " --xmloutput " + fileInfoToPath info 
            | IncludeDecoysInXML         -> " --decoy-xml-output " 
        
        type SVMFeatureOptions =
            /// Output final SVM weights to given file. (one per line)
            | OUT_SVMWeights of System.IO.FileInfo
            /// Read initial SVM weights from given file (one per line)
            | IN_SVMWeights of System.IO.FileInfo
            /// Use given feature name as initial search direction, can be negated to indicate that a lower value is better.
            //TODO:
            //| FEATURENAME of System.IO.FileInfo
            /// Use unit normalization [0-1] on features instead of standard deviation normalization.
            | UnitNorm
            /// Override error check and do not fall back on default score vector in case of suspect score vector from SVM.
            | Override
            /// nclude description of correct features, i.e. features describing the difference between the observed and 
            /// predicted isoelectric point, retention time and precursor mass. See this page for a more detailed description
            | DOC
            /// Retention time features are calculated as in Klammer et al. instead of with Elude. Only available if -D is set.
            | Klammer
        
        let private stringOfSVMFO fileInfoToPath (sVMFO: SVMFeatureOptions) = 
            match sVMFO with 
            | OUT_SVMWeights info   -> " --weights " + fileInfoToPath info 
            | IN_SVMWeights info    -> " --init-weights " + fileInfoToPath info 
            // TODO
            // | FEATURENAME info          -> "--results-psms " + info.FullName
            | UnitNorm              -> " --unitnorm "
            | Override              -> " --override "
            | DOC                   -> " --doc " 
            | Klammer               -> " --klammer "  

        type SVMTrainingOptions =
            /// Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other PSMs. 
            /// Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal. Default = 0.
            | SubsetTraining of float
            /// Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.
            | Cpos of float
            /// Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified or if -p is not specified
            | Cneg of float
            /// False discovery rate threshold for evaluating best cross validation result and reported end result. Default = 0.01.
            | FDR_CrossValidation of float 
            /// False discovery rate threshold to define positive examples in training. Set to testFDR if 0. Default = 0.01.
            | FDR_PositiveExamples of float 
            /// Maximal number of iterations
            | MaxIterations of int
            /// Quicker execution by reduced internal cross-validation.
            | QuickValidation
            /// Report performance on test set each iteration.
            | ReportPerfomanceAfterIteration
            /// Set seed of the random number generator. Default = 1
            | SeedRndNumberGenerator of float 

        let private stringOfSVMTO (sVMTO: SVMTrainingOptions) = 
            match sVMTO with 
            | SubsetTraining v               -> " --subset-max-train " + string v
            | Cpos v                         -> " --Cpos "  + string v
            | Cneg v                         -> " --Cneg" + string v
            | FDR_CrossValidation v          -> " --testFDR " + string v
            | FDR_PositiveExamples v         -> " --trainFDR " + string v
            | MaxIterations n                 -> " --maxiter " + string n 
            | QuickValidation                -> " --quick-validation "
            | ReportPerfomanceAfterIteration -> " --test-each-iteration "
            | SeedRndNumberGenerator v       -> " --seed" + string v

        type ProteinInferenceOptions_Percolator = 
            /// Use the picked protein-level FDR to infer protein probabilities. Provide the fasta file as the argument to this flag, 
            /// which will be used for protein grouping based on an in-silico digest. If no fasta file is available or protein grouping 
            /// is not desired, set this flag to "auto" to skip protein grouping.
            | Fasta of System.IO.FileInfo
            /// Define the text pattern to identify decoy proteins in the database. Default = "random_".
            | ProteinDecoyPattern of string 
            /// Type of enzyme "no_enzyme","elastase","pepsin","proteinasek","thermolysin","trypsinp","chymotrypsin","lys-n","lys-c",
            /// "arg-c","asp-n","glu-c","trypsin" default="trypsin"
            | Protease of string 
            /// By default, if the peptides associated with protein A are a proper subset of the peptides associated with protein B, 
            /// then protein A is eliminated and all the peptides are considered as evidence for protein B. Note that this filtering is 
            /// done based on the complete set of peptides in the database, not based on the identified peptides in the search results.
            /// Alternatively, if this option is set and if all of the identified peptides associated with protein B are also associated 
            /// with protein A, then Percolator will report a comma-separated list of protein IDs, where the full-length protein B is first
            /// in the list and the fragment protein A is listed second. Not available for Fido.
            | ReportProteinFragments
            /// If multiple database proteins contain exactly the same set of peptides, then Percolator will randomly discard all but one of the proteins.
            /// If this option is set, then the IDs of these duplicated proteins will be reported as a comma-separated list.
            /// Not available for Fido.
            | ReportProteinDuplicates

        let private stringOfPI fileInfoToPath (pI:ProteinInferenceOptions_Percolator) =
            match pI with 
            | Fasta info                    -> " --picked-protein " + fileInfoToPath info 
            | ProteinDecoyPattern pattern   -> " --protein-decoy-pattern " + pattern
            | Protease protease             -> " --protein-enzyme " + protease
            | ReportProteinFragments        -> " --protein-report-fragments " 
            | ReportProteinDuplicates       -> " --protein-report-duplicates "         

        type ProteinInferenceOptions_FIDO = 
            /// Use the Fido algorithm to infer protein probabilities.
            | UseFido
            /// Set Fido's probability with which a present protein emits an associated peptide. Set by grid search if not specified.
            | Alpha of float 
            /// Set Fido's probability of creation of a peptide from noise. Set by grid search if not specified.
            | Beta of float
            /// Set Fido's prior probability that a protein is present in the sample. Set by grid search if not specified.
            | Gamma of float 
            /// Estimate empirical p-values and q-values using target-decoy analysis.
            | EmpricialQValue
            /// Q-value threshold that will be used in the computation of the MSE and ROC AUC score in the grid search.
            /// Recommended 0.05 for normal size datasets and 0.1 for big size datasets. Default = 0.1.
            | QValueThreshold of float 
    
        //type SpeedUpOptions_FIDO =     
            /// Setting the gridsearch-depth to 0 (fastest), 1 or 2 (slowest) controls how much computational time is required
            /// for the estimation of alpha, beta and gamma parameters for Fido. Default = 0.
            | GridSearchDepth of float 
            /// Apply the specified threshold to PSM, peptide and protein probabilities to obtain a faster estimate of the alpha, 
            /// beta and gamma parameters. Default = 0; Recommended when set = 0.2.
            | GridSearchSpeed of float 
            /// Do not approximate the posterior distribution by allowing large graph components to be split into subgraphs. 
            /// The splitting is done by duplicating peptides with low probabilities. Splitting continues until the number of 
            /// possible configurations of each subgraph is below 2^18.
            | NoSubgraphSplitting
            /// To speed up inference, proteins for which none of the associated peptides has a probability exceeding 
            /// the specified threshold will be assigned probability = 0. Default = 0.01.
            | ProteinTruncationThreshold 

        let private stringOfPIF (pIF:ProteinInferenceOptions_FIDO) =
            match pIF with 
            | UseFido                -> " --fido-protein "
            | Alpha value            -> " --fido-alpha " + string value
            | Beta value             -> " --fido-beta " + string value
            | Gamma value            -> " --fido-gamma " + string value
            | EmpricialQValue        -> " --fido-empirical-protein-q " 
            | QValueThreshold value  -> " --fido-gridsearch-mse-threshold " + string value 

        //let private stringOfSUF (sUF:SpeedUpOptions_FIDO) =
        //    match sUF with 
            | GridSearchDepth value      -> " --fido-gridsearch-depth " + string value
            | GridSearchSpeed value      -> " --fido-fast-gridsearch " + string value
            | NoSubgraphSplitting        -> " --fido-no-split-large-components " 
            | ProteinTruncationThreshold -> " --fido-protein-truncation-threshold " 
            

        type PercolatorParams =
            | GeneralOptions            of seq<GeneralOptions>
            | FileInputOptions          of seq<FileInputOptions>
            | FileOutputOptions         of seq<FileOutputOptions> 
            | SVMFeatureOptions         of seq<SVMFeatureOptions>
            | SVMTrainingOptions        of seq<SVMTrainingOptions>
            | ProteinInferenceOptions_FIDO  of seq<ProteinInferenceOptions_FIDO>
            | ProteinInferenceOptions_Percolator of seq<ProteinInferenceOptions_Percolator>
    
        let stringOf fileInfoToPath (p:PercolatorParams) = 
            let iterCustom f s =
                Seq.map f s
                |> String.concat ""
            match p with
            | GeneralOptions                     s -> iterCustom stringOfgO s
            | FileInputOptions                   s -> iterCustom (stringOfFIO fileInfoToPath) s
            | FileOutputOptions                  s -> iterCustom (stringOfFOO fileInfoToPath) s
            | SVMFeatureOptions                  s -> iterCustom (stringOfSVMFO fileInfoToPath) s
            | SVMTrainingOptions                 s -> iterCustom stringOfSVMTO s
            | ProteinInferenceOptions_FIDO       s -> iterCustom stringOfPIF s  
            | ProteinInferenceOptions_Percolator s -> iterCustom (stringOfPI fileInfoToPath) s 

    open System 
    open System.Diagnostics

    type OperatingSystem = 
        | Windows
        | Ubuntu 

    type PercolatorWrapper(os:OperatingSystem,?exePath) = 
        let createArguments (f : 'a -> string) (ps:seq<'a>) =
            ps |> Seq.map f
            |> String.concat " "

        let createUnixShellProcess name (f : 'TParam -> string) (arg:string) =
            let path = 
                match exePath with 
                | Some path -> path 
                | None      -> """-c "/usr/bin/percolator """
            let beginTime = DateTime.UtcNow
            printfn "Starting %s..." name
            let p =                        
                new ProcessStartInfo
                 (FileName = "bash.exe", UseShellExecute = false, Arguments = path + arg + "\"",
                  RedirectStandardError = false, CreateNoWindow = false,
                  RedirectStandardOutput = false, RedirectStandardInput = false)
                |> Process.Start
            p.WaitForExit()
            p.Close()
            printfn "%s done." name
            printfn "Elapsed time: %A" (beginTime.Subtract(DateTime.UtcNow))

        let createWindowsProcess name (f : 'TParam -> string) (arg:string) =           
            let exePath' = 
                match exePath with 
                | Some path -> path 
                | None      -> """-c "/usr/bin/percolator """
            let beginTime = DateTime.UtcNow
            printfn "Starting %s..." name
            let p =                        
                new ProcessStartInfo
                  (FileName = exePath', UseShellExecute = false, Arguments = arg, 
                   RedirectStandardError = false, CreateNoWindow = true, 
                   RedirectStandardOutput = false, RedirectStandardInput = true) 
                |> Process.Start
            p.WaitForExit()
            p.Close()
            printfn "%s done." name
            printfn "Elapsed time: %A" (beginTime.Subtract(DateTime.UtcNow))
        
            
        
        member this.Percolate (ps:seq<Parameters.PercolatorParams>) =
            match os with 
            | Windows ->
                let arg = createArguments (Parameters.stringOf Parameters.fileInfoToWindowsPath) ps 
                createWindowsProcess (sprintf "Percolator with %s" arg) (Parameters.stringOf Parameters.fileInfoToWindowsPath) arg

            | Ubuntu -> 
                let arg = createArguments (Parameters.stringOf Parameters.fileInfoToLinuxPath) ps 
                createUnixShellProcess (sprintf "Percolator with %s" arg) (Parameters.stringOf Parameters.fileInfoToLinuxPath) arg

      