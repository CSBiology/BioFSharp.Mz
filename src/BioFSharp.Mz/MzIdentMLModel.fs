namespace BioFSharp.Mz


open BioFSharp
open BioFSharp.IO

open System
open FSharp.Care
open FSharp.Care.Collections
open FSharp.Care.Monads
open AminoAcids 
open ModificationInfo
open System.Data
open System.Data.SQLite


/// Module containing all data types used in the model
module DataModel =

    /// Term (cv) accession number
    type TermId = string

    /// Define a controlled vocabulary term
    [<CustomEquality; CustomComparison>]
    type Term = {
        /// accession number
        Id          : TermId
        FK_Ontology : string
        Name        : string
        RowVersion  : System.DateTime
        }
        with
            override this.Equals(other) =
                match other with
                    | :? Term as o -> this.Id = o.Id && this.RowVersion = o.RowVersion
                    | _ -> false

                override this.GetHashCode() = hash this.Id
                            
                interface System.IComparable with
                    member this.CompareTo other =
                        match other with
                        | :? Term as o -> compare this.Id o.Id
                        | _ -> invalidArg "Other" "cannot compare values of different types"




    /// Define a controlled vocabulary parameter
    [<CustomEquality; CustomComparison>]
    type CvParam = {
        Id          : System.Guid
        // FK_ParamContainer   : System.Guid
        // FK_Term             : System.Guid
        Term        : Term
        Value       : System.IConvertible
        Unit        : Term option
        RowVersion  : System.DateTime
        }    
        with
            override this.Equals(other) =
                match other with
                    | :? CvParam as o -> this.Id = o.Id && this.RowVersion = o.RowVersion
                    | _ -> false

                override this.GetHashCode() = hash this.Id
                            
                interface System.IComparable with
                    member this.CompareTo other =
                        match other with
                        | :? CvParam as o -> compare this.Term.Id o.Term.Id
                        | _ -> invalidArg "Other" "cannot compare values of different types"

/// Module containing functions to operate on controlled vocabulary term
module Term =

    open DataModel

    /// Creates a controlled vocabulary term with rowVersion
    let initOf id ontologyId name rowVersion =
        {Id=id;FK_Ontology=ontologyId;Name=name;RowVersion=rowVersion}

    /// Creates a controlled vocabulary term
    let create id ontologyId name =
        {Id=id;FK_Ontology=ontologyId;Name=name;RowVersion=System.DateTime.Now}


/// Module containing functions to operate on controlled vocabulary parameter
module CvParam =

    open DataModel
    

    /// Creates a controlled vocabulary parameter
    let create id term value =
        {Id=id;Term=term;Value=value;Unit=None;RowVersion=System.DateTime.Now}

    /// Creates a controlled vocabulary parameter with unit term 
    let createWithUnit id term value unit =
        {Id=id;Term=term;Value=value;Unit=Some unit;RowVersion=System.DateTime.Now}

    /// Creates a controlled vocabulary parameter value
    let ofValue term value =
        {Id=System.Guid.NewGuid();Term=term;Value=value;Unit=None;RowVersion=System.DateTime.Now}

    /// Creates a controlled vocabulary parameter value with unit term 
    let createCvParamValueWithUnit term value unit =
        {Id=System.Guid.NewGuid();Term=term;Value=value;Unit=Some unit;RowVersion=System.DateTime.Now}

//    let ofValue termId ontologyId name value =
//        let term = Term.create termId ontologyId name
//        ofValue term value


/// Module containing functions to operate on ParamContainer
module ParamContainer =
    open DataModel

    /// Type ParamContainer
    type ParamContainer = System.Collections.Generic.Dictionary<TermId,CvParam>

    let addOrUpdateInPlace (param:CvParam) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(param.Term.Id) then
            paramContainer.[param.Term.Id] <- param 
        else
            paramContainer. Add(param.Term.Id,param)
        paramContainer

    let addOrUpdateInPlaceOf (cvParams:seq<CvParam>) (paramContainer:ParamContainer) =
        cvParams |> Seq.iter (fun p -> addOrUpdateInPlace p paramContainer |> ignore)
        paramContainer

    let ofSeq (cvParams:seq<CvParam>)  =
        let dict = System.Collections.Generic.Dictionary<TermId,CvParam>()
        cvParams |> Seq.iter (fun p -> addOrUpdateInPlace p dict |> ignore) 
        dict

    /// Get CvParam by id (Function can fail)
    let getCvParam (paramId:TermId) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(paramId) then
           paramContainer.[paramId]
        else
            failwithf "Parameter not found %s" paramId
    
    /// Get optional CvParam by id 
    let tryGetCvParam (paramId:TermId) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(paramId) then
            Some paramContainer.[paramId]
        else
            None

    /// Get value as float by termid or 'nan'
    let getValueAsFloat (paramId:TermId) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(paramId) then
            System.Convert.ToDouble (paramContainer.[paramId].Value)
        else
            nan

    /// Get value as int32 by termid or '-1'
    let getValueAsInt (paramId:TermId) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(paramId) then
            System.Convert.ToInt32 (paramContainer.[paramId].Value)
        else
            -1

    /// Get value as string by termid or empty string
    let getValueAsString (paramId:TermId) (paramContainer:ParamContainer) =
        if paramContainer.ContainsKey(paramId) then
            System.Convert.ToString (paramContainer.[paramId].Value)
        else
            ""
module Db =

    module DBSequence =
        open DataModel
        open ParamContainer


        type DBSequence = {
            ID : int 
            Accession : string 
            Name : string 
            SearchDBID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createDBSequence id accession name searchdbid rowversion  = {
            ID=id; Accession=accession; Name=name; SearchDBID=searchdbid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createDBSequenceWith id accession name searchdbid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Accession=accession; Name=name; SearchDBID=searchdbid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (dbsequence:DBSequence)  = 
            ParamContainer.addOrUpdateInPlace param dbsequence.ParamContainer |> ignore
            dbsequence

        let createDBSequenceTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE DBSequence (
                ID INTEGER NOT NULL,
                Accession TEXT NOT NULL,
                Name TEXT NOT NULL,
                SearchDBID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createDBSequenceParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE DBSequenceParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertDBSequence (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO DBSequence (
                    ID,
                    Accession,
                    Name,
                    SearchDBID,
                    RowVersion)
                    VALUES (
                        @id,
                        @accession,
                        @name,
                        @searchdbid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@accession",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@searchdbid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id accession name searchdbid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@accession"].Value <- accession
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@searchdbid"].Value <- searchdbid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertDBSequenceParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO DBSequenceParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertDBSequenceToDb insertInTableF insertInParamTableF dbsequence =
             insertInTableF
                dbsequence.ID
                dbsequence.Accession
                dbsequence.Name
                dbsequence.SearchDBID
                dbsequence.RowVersion
             insertInParamTableF dbsequence.ParamContainer

            ///
        let prepareSelectDBSequencebyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMDBSequence WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectDBSequencebyAccession (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMDBSequence WHERE Accession=@accession" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ accession ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun accession->
            cmd.Parameters.["@ accession "].Value <- accession
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectDBSequencebyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMDBSequence WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectDBSequencebySearchDBID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMDBSequence WHERE SearchDBID=@searchdbid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ searchdbid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun searchdbid->
            cmd.Parameters.["@ searchdbid "].Value <- searchdbid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectDBSequencebyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMDBSequence WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module PeptideEvidence =
        open DataModel
        open ParamContainer


        type PeptideEvidence = {
            ID : int 
            DBSequenceID : int 
            PeptideID : int 
            isDecoy : string option;
            Frame : string option;
            Start : int option;
            End : int option;
            Pre : string option;
            Post : string option;
            TranslationsID : int option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createPeptideEvidence id dbsequenceid peptideid isdecoy frame start endOf pre post translationsid rowversion  = {
            ID=id; DBSequenceID=dbsequenceid; PeptideID=peptideid; isDecoy=isdecoy; Frame=frame; Start=start; End=endOf; Pre=pre; Post=post; TranslationsID=translationsid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createPeptideEvidenceWith id dbsequenceid peptideid isdecoy frame start endOf pre post translationsid rowversion (cvParams:seq<CvParam>) = {
            ID=id; DBSequenceID=dbsequenceid; PeptideID=peptideid; isDecoy=isdecoy; Frame=frame; Start=start; End=endOf; Pre=pre; Post=post; TranslationsID=translationsid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (peptideevidence:PeptideEvidence)  = 
            ParamContainer.addOrUpdateInPlace param peptideevidence.ParamContainer |> ignore
            peptideevidence

        let createPeptideEvidenceTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PeptideEvidence (
                ID INTEGER NOT NULL,
                DBSequenceID INTEGER NOT NULL,
                PeptideID INTEGER NOT NULL,
                isDecoy TEXT,
                Frame TEXT,
                Start INTEGER,
                End INTEGER,
                Pre TEXT,
                Post TEXT,
                TranslationsID INTEGER,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT fk_DBSequence_id FOREIGN KEY (DBSequenceID) REFERENCES DBSequence (ID),
                CONSTRAINT fk_Peptide_id FOREIGN KEY (PeptideID) REFERENCES Peptide (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createPeptideEvidenceParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PeptideEvidenceParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertPeptideEvidence (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PeptideEvidence (
                    ID,
                    DBSequenceID,
                    PeptideID,
                    isDecoy,
                    Frame,
                    Start,
                    End,
                    Pre,
                    Post,
                    TranslationsID,
                    RowVersion)
                    VALUES (
                        @id,
                        @dbsequenceid,
                        @peptideid,
                        @isdecoy,
                        @frame,
                        @start,
                        @end,
                        @pre,
                        @post,
                        @translationsid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@dbsequenceid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@peptideid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@isdecoy",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@frame",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@start",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@end",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@pre",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@post",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@translationsid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id dbsequenceid peptideid isdecoy frame start endOf pre post translationsid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@dbsequenceid"].Value <- dbsequenceid
                cmd.Parameters.["@peptideid"].Value <- peptideid
                cmd.Parameters.["@isdecoy"].Value <- isdecoy
                cmd.Parameters.["@frame"].Value <- frame
                cmd.Parameters.["@start"].Value <- start
                cmd.Parameters.["@end"].Value <- endOf
                cmd.Parameters.["@pre"].Value <- pre
                cmd.Parameters.["@post"].Value <- post
                cmd.Parameters.["@translationsid"].Value <- translationsid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertPeptideEvidenceParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PeptideEvidenceParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertPeptideEvidenceToDb insertInTableF insertInParamTableF peptideevidence =
             insertInTableF
                peptideevidence.ID
                peptideevidence.DBSequenceID
                peptideevidence.PeptideID
                peptideevidence.isDecoy
                peptideevidence.Frame
                peptideevidence.Start
                peptideevidence.End
                peptideevidence.Pre
                peptideevidence.Post
                peptideevidence.TranslationsID
                peptideevidence.RowVersion
             insertInParamTableF peptideevidence.ParamContainer

            ///
        let prepareSelectPeptideEvidencebyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyDBSequenceID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE DBSequenceID=@dbsequenceid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ dbsequenceid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun dbsequenceid->
            cmd.Parameters.["@ dbsequenceid "].Value <- dbsequenceid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyPeptideID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE PeptideID=@peptideid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ peptideid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun peptideid->
            cmd.Parameters.["@ peptideid "].Value <- peptideid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyisDecoy (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE isDecoy=@isdecoy" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ isdecoy ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun isdecoy->
            cmd.Parameters.["@ isdecoy "].Value <- isdecoy
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyFrame (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE Frame=@frame" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ frame ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun frame->
            cmd.Parameters.["@ frame "].Value <- frame
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyStart (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE Start=@start" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ start ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun start->
            cmd.Parameters.["@ start "].Value <- start
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyEnd (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE End=@end" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ end ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun endOf ->
            cmd.Parameters.["@ end "].Value <- endOf
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyPre (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE Pre=@pre" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ pre ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun pre->
            cmd.Parameters.["@ pre "].Value <- pre
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyPost (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE Post=@post" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ post ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun post->
            cmd.Parameters.["@ post "].Value <- post
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyTranslationsID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE TranslationsID=@translationsid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ translationsid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun translationsid->
            cmd.Parameters.["@ translationsid "].Value <- translationsid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideEvidencebyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideEvidence WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetInt32(5), reader.GetInt32(6), reader.GetString(7), reader.GetString(8), reader.GetInt32(9), reader.GetDateTime(10) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Peptide =
        open DataModel
        open ParamContainer


        type Peptide = {
            ID : string 
            Sequence : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createPeptide id sequence rowversion  = {
            ID=id; Sequence=sequence; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createPeptideWith id sequence rowversion (cvParams:seq<CvParam>) = {
            ID=id; Sequence=sequence; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (peptide:Peptide)  = 
            ParamContainer.addOrUpdateInPlace param peptide.ParamContainer |> ignore
            peptide

        let createPeptideTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Peptide (
                ID TEXT NOT NULL,
                Sequence TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createPeptideParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PeptideParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertPeptide (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Peptide (
                    ID,
                    Sequence,
                    RowVersion)
                    VALUES (
                        @id,
                        @sequence,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@sequence",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id sequence rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@sequence"].Value <- sequence
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertPeptideParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PeptideParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertPeptideToDb insertInTableF insertInParamTableF peptide =
             insertInTableF
                peptide.ID
                peptide.Sequence
                peptide.RowVersion
             insertInParamTableF peptide.ParamContainer

            ///
        let prepareSelectPeptidebyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptide WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetString(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptidebySequence (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptide WHERE Sequence=@sequence" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ sequence ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetString(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun sequence->
            cmd.Parameters.["@ sequence "].Value <- sequence
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptidebyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptide WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetString(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module SpectrumIdentificationItem =
        open DataModel
        open ParamContainer


        type SpectrumIdentificationItem = {
            ID : int 
            SpectrumIdentificationResultID : int option;
            SampleID : int option;
            PeptideID : int 
            MassTableID : int option;
            Name : string option;
            PassThreshold : string 
            Rank : int option;
            CalculatedMassToCharge : float option;
            ExperimentalMassToCharge : float 
            ChargeState : int 
            CalculatedPI : float option;
            Fragmentation : System.DateTime option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createSpectrumIdentificationItem id spectrumidentificationresultid sampleid peptideid masstableid name passthreshold rank calculatedmasstocharge experimentalmasstocharge chargestate calculatedpi fragmentation rowversion  = {
            ID=id; SpectrumIdentificationResultID=spectrumidentificationresultid; SampleID=sampleid; PeptideID=peptideid; MassTableID=masstableid; Name=name; PassThreshold=passthreshold; Rank=rank; CalculatedMassToCharge=calculatedmasstocharge; ExperimentalMassToCharge=experimentalmasstocharge; ChargeState=chargestate; CalculatedPI=calculatedpi; Fragmentation=fragmentation; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createSpectrumIdentificationItemWith id spectrumidentificationresultid sampleid peptideid masstableid name passthreshold rank calculatedmasstocharge experimentalmasstocharge chargestate calculatedpi fragmentation rowversion (cvParams:seq<CvParam>) = {
            ID=id; SpectrumIdentificationResultID=spectrumidentificationresultid; SampleID=sampleid; PeptideID=peptideid; MassTableID=masstableid; Name=name; PassThreshold=passthreshold; Rank=rank; CalculatedMassToCharge=calculatedmasstocharge; ExperimentalMassToCharge=experimentalmasstocharge; ChargeState=chargestate; CalculatedPI=calculatedpi; Fragmentation=fragmentation; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (spectrumidentificationitem:SpectrumIdentificationItem)  = 
            ParamContainer.addOrUpdateInPlace param spectrumidentificationitem.ParamContainer |> ignore
            spectrumidentificationitem

        let createSpectrumIdentificationItemTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationItem (
                ID INTEGER NOT NULL,
                SpectrumIdentificationResultID INTEGER,
                SampleID INTEGER,
                PeptideID INTEGER NOT NULL,
                MassTableID INTEGER,
                Name TEXT,
                PassThreshold TEXT NOT NULL,
                Rank INTEGER,
                CalculatedMassToCharge REAL,
                ExperimentalMassToCharge REAL NOT NULL,
                ChargeState INTEGER NOT NULL,
                CalculatedPI REAL,
                Fragmentation BLOB,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_SpecIdItem_Peptide_ID FOREIGN KEY (PeptideID) REFERENCES Peptide (ID),
                CONSTRAINT FK_SpecIDItem_SpecIDRes_ID FOREIGN KEY (SpectrumIdentificationResultID) REFERENCES SpectrumIdentificationResult (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createSpectrumIdentificationItemParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationItemParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertSpectrumIdentificationItem (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationItem (
                    ID,
                    SpectrumIdentificationResultID,
                    SampleID,
                    PeptideID,
                    MassTableID,
                    Name,
                    PassThreshold,
                    Rank,
                    CalculatedMassToCharge,
                    ExperimentalMassToCharge,
                    ChargeState,
                    CalculatedPI,
                    Fragmentation,
                    RowVersion)
                    VALUES (
                        @id,
                        @spectrumidentificationresultid,
                        @sampleid,
                        @peptideid,
                        @masstableid,
                        @name,
                        @passthreshold,
                        @rank,
                        @calculatedmasstocharge,
                        @experimentalmasstocharge,
                        @chargestate,
                        @calculatedpi,
                        @fragmentation,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@spectrumidentificationresultid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@sampleid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@peptideid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@masstableid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@passthreshold",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rank",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@calculatedmasstocharge",Data.DbType.Double) |> ignore
            cmd.Parameters.Add("@experimentalmasstocharge",Data.DbType.Double) |> ignore
            cmd.Parameters.Add("@chargestate",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@calculatedpi",Data.DbType.Double) |> ignore
            cmd.Parameters.Add("@fragmentation",Data.DbType.DateTime) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id spectrumidentificationresultid sampleid peptideid masstableid name passthreshold rank calculatedmasstocharge experimentalmasstocharge chargestate calculatedpi fragmentation rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@spectrumidentificationresultid"].Value <- spectrumidentificationresultid
                cmd.Parameters.["@sampleid"].Value <- sampleid
                cmd.Parameters.["@peptideid"].Value <- peptideid
                cmd.Parameters.["@masstableid"].Value <- masstableid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@passthreshold"].Value <- passthreshold
                cmd.Parameters.["@rank"].Value <- rank
                cmd.Parameters.["@calculatedmasstocharge"].Value <- calculatedmasstocharge
                cmd.Parameters.["@experimentalmasstocharge"].Value <- experimentalmasstocharge
                cmd.Parameters.["@chargestate"].Value <- chargestate
                cmd.Parameters.["@calculatedpi"].Value <- calculatedpi
                cmd.Parameters.["@fragmentation"].Value <- fragmentation
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertSpectrumIdentificationItemParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationItemParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertSpectrumIdentificationItemToDb insertInTableF insertInParamTableF spectrumidentificationitem =
             insertInTableF
                spectrumidentificationitem.ID
                spectrumidentificationitem.SpectrumIdentificationResultID
                spectrumidentificationitem.SampleID
                spectrumidentificationitem.PeptideID
                spectrumidentificationitem.MassTableID
                spectrumidentificationitem.Name
                spectrumidentificationitem.PassThreshold
                spectrumidentificationitem.Rank
                spectrumidentificationitem.CalculatedMassToCharge
                spectrumidentificationitem.ExperimentalMassToCharge
                spectrumidentificationitem.ChargeState
                spectrumidentificationitem.CalculatedPI
                spectrumidentificationitem.Fragmentation
                spectrumidentificationitem.RowVersion
             insertInParamTableF spectrumidentificationitem.ParamContainer

            ///
        let prepareSelectSpectrumIdentificationItembyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembySpectrumIdentificationResultID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE SpectrumIdentificationResultID=@spectrumidentificationresultid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectrumidentificationresultid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun spectrumidentificationresultid->
            cmd.Parameters.["@ spectrumidentificationresultid "].Value <- spectrumidentificationresultid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembySampleID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE SampleID=@sampleid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ sampleid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun sampleid->
            cmd.Parameters.["@ sampleid "].Value <- sampleid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyPeptideID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE PeptideID=@peptideid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ peptideid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun peptideid->
            cmd.Parameters.["@ peptideid "].Value <- peptideid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyMassTableID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE MassTableID=@masstableid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ masstableid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun masstableid->
            cmd.Parameters.["@ masstableid "].Value <- masstableid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyPassThreshold (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE PassThreshold=@passthreshold" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ passthreshold ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun passthreshold->
            cmd.Parameters.["@ passthreshold "].Value <- passthreshold
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyRank (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE Rank=@rank" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rank ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun rank->
            cmd.Parameters.["@ rank "].Value <- rank
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyCalculatedMassToCharge (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE CalculatedMassToCharge=@calculatedmasstocharge" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ calculatedmasstocharge ",Data.DbType.Double) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun calculatedmasstocharge->
            cmd.Parameters.["@ calculatedmasstocharge "].Value <- calculatedmasstocharge
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyExperimentalMassToCharge (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE ExperimentalMassToCharge=@experimentalmasstocharge" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ experimentalmasstocharge ",Data.DbType.Double) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun experimentalmasstocharge->
            cmd.Parameters.["@ experimentalmasstocharge "].Value <- experimentalmasstocharge
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyChargeState (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE ChargeState=@chargestate" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ chargestate ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun chargestate->
            cmd.Parameters.["@ chargestate "].Value <- chargestate
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyCalculatedPI (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE CalculatedPI=@calculatedpi" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ calculatedpi ",Data.DbType.Double) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun calculatedpi->
            cmd.Parameters.["@ calculatedpi "].Value <- calculatedpi
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyFragmentation (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE Fragmentation=@fragmentation" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ fragmentation ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun fragmentation->
            cmd.Parameters.["@ fragmentation "].Value <- fragmentation
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationItembyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationItem WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetString(5), reader.GetString(6), reader.GetInt32(7), reader.GetDouble(8), reader.GetDouble(9), reader.GetInt32(10), reader.GetDouble(11), reader.GetDateTime(12), reader.GetDateTime(13) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Term =
        open DataModel
        open ParamContainer


        type Term = {
            ID : int 
            OntologyID : int 
            Name : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createTerm id ontologyid name rowversion  = {
            ID=id; OntologyID=ontologyid; Name=name; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createTermWith id ontologyid name rowversion (cvParams:seq<CvParam>) = {
            ID=id; OntologyID=ontologyid; Name=name; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (term:Term)  = 
            ParamContainer.addOrUpdateInPlace param term.ParamContainer |> ignore
            term

        let createTermTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Term (
                ID INTEGER NOT NULL,
                OntologyID INTEGER NOT NULL,
                Name TEXT NOT NULL,
                RowVersion  BLOB(8) NOT NULL DEFAULT 0,
                CONSTRAINT FK_Term_Ontology_ID FOREIGN KEY (OntologyID) REFERENCES Ontology (Name)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createTermParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE TermParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertTerm (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Term (
                    ID,
                    OntologyID,
                    Name,
                    RowVersion)
                    VALUES (
                        @id,
                        @ontologyid,
                        @name,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@ontologyid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id ontologyid name rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@ontologyid"].Value <- ontologyid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertTermParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO TermParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertTermToDb insertInTableF insertInParamTableF term =
             insertInTableF
                term.ID
                term.OntologyID
                term.Name
                term.RowVersion
             insertInParamTableF term.ParamContainer

            ///
        let prepareSelectTermbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTerm WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermbyOntologyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTerm WHERE OntologyID=@ontologyid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ ontologyid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun ontologyid->
            cmd.Parameters.["@ ontologyid "].Value <- ontologyid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTerm WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTerm WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module SpectrumIdentificationResult =
        open DataModel
        open ParamContainer


        type SpectrumIdentificationResult = {
            ID : int 
            SpectrumID : string 
            SpectraDataID : string 
            SpectrumIdentficationListID : int option;
            Name : string option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createSpectrumIdentificationResult id spectrumid spectradataid spectrumidentficationlistid name rowversion  = {
            ID=id; SpectrumID=spectrumid; SpectraDataID=spectradataid; SpectrumIdentficationListID=spectrumidentficationlistid; Name=name; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createSpectrumIdentificationResultWith id spectrumid spectradataid spectrumidentficationlistid name rowversion (cvParams:seq<CvParam>) = {
            ID=id; SpectrumID=spectrumid; SpectraDataID=spectradataid; SpectrumIdentficationListID=spectrumidentficationlistid; Name=name; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (spectrumidentificationresult:SpectrumIdentificationResult)  = 
            ParamContainer.addOrUpdateInPlace param spectrumidentificationresult.ParamContainer |> ignore
            spectrumidentificationresult

        let createSpectrumIdentificationResultTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationResult (
                ID INTEGER NOT NULL,
                SpectrumID TEXT NOT NULL,
                SpectraDataID TEXT NOT NULL,
                SpectrumIdentficationListID INTEGER,
                Name TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_SpecIdentResults_SpecIdentList_ID FOREIGN KEY (SpectrumIdentficationListID) REFERENCES SpectrumIdentificationList (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createSpectrumIdentificationResultParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationResultParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertSpectrumIdentificationResult (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationResult (
                    ID,
                    SpectrumID,
                    SpectraDataID,
                    SpectrumIdentficationListID,
                    Name,
                    RowVersion)
                    VALUES (
                        @id,
                        @spectrumid,
                        @spectradataid,
                        @spectrumidentficationlistid,
                        @name,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@spectrumid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@spectradataid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@spectrumidentficationlistid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id spectrumid spectradataid spectrumidentficationlistid name rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@spectrumid"].Value <- spectrumid
                cmd.Parameters.["@spectradataid"].Value <- spectradataid
                cmd.Parameters.["@spectrumidentficationlistid"].Value <- spectrumidentficationlistid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertSpectrumIdentificationResultParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationResultParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertSpectrumIdentificationResultToDb insertInTableF insertInParamTableF spectrumidentificationresult =
             insertInTableF
                spectrumidentificationresult.ID
                spectrumidentificationresult.SpectrumID
                spectrumidentificationresult.SpectraDataID
                spectrumidentificationresult.SpectrumIdentficationListID
                spectrumidentificationresult.Name
                spectrumidentificationresult.RowVersion
             insertInParamTableF spectrumidentificationresult.ParamContainer

            ///
        let prepareSelectSpectrumIdentificationResultbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationResultbySpectrumID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE SpectrumID=@spectrumid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectrumid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun spectrumid->
            cmd.Parameters.["@ spectrumid "].Value <- spectrumid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationResultbySpectraDataID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE SpectraDataID=@spectradataid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectradataid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun spectradataid->
            cmd.Parameters.["@ spectradataid "].Value <- spectradataid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationResultbySpectrumIdentficationListID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE SpectrumIdentficationListID=@spectrumidentficationlistid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectrumidentficationlistid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun spectrumidentficationlistid->
            cmd.Parameters.["@ spectrumidentficationlistid "].Value <- spectrumidentficationlistid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationResultbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationResultbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationResult WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module SpectrumIdentificationList =
        open DataModel
        open ParamContainer


        type SpectrumIdentificationList = {
            ID : int 
            Name : string option;
            NumSequencesSearched : int option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createSpectrumIdentificationList id name numsequencessearched rowversion  = {
            ID=id; Name=name; NumSequencesSearched=numsequencessearched; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createSpectrumIdentificationListWith id name numsequencessearched rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; NumSequencesSearched=numsequencessearched; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (spectrumidentificationlist:SpectrumIdentificationList)  = 
            ParamContainer.addOrUpdateInPlace param spectrumidentificationlist.ParamContainer |> ignore
            spectrumidentificationlist

        let createSpectrumIdentificationListTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationList (
                ID INTEGER NOT NULL,
                Name TEXT,
                NumSequencesSearched INTEGER,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createSpectrumIdentificationListParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationListParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertSpectrumIdentificationList (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationList (
                    ID,
                    Name,
                    NumSequencesSearched,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @numsequencessearched,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@numsequencessearched",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name numsequencessearched rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@numsequencessearched"].Value <- numsequencessearched
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertSpectrumIdentificationListParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationListParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertSpectrumIdentificationListToDb insertInTableF insertInParamTableF spectrumidentificationlist =
             insertInTableF
                spectrumidentificationlist.ID
                spectrumidentificationlist.Name
                spectrumidentificationlist.NumSequencesSearched
                spectrumidentificationlist.RowVersion
             insertInParamTableF spectrumidentificationlist.ParamContainer

            ///
        let prepareSelectSpectrumIdentificationListbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationList WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationListbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationList WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationListbyNumSequencesSearched (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationList WHERE NumSequencesSearched=@numsequencessearched" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ numsequencessearched ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun numsequencessearched->
            cmd.Parameters.["@ numsequencessearched "].Value <- numsequencessearched
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationListbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationList WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ModLocation =
        open DataModel
        open ParamContainer


        type ModLocation = {
            ID : int 
            PeptideID : int 
            ModificationID : int 
            Location : int 
            Residue : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createModLocation id peptideid modificationid location residue rowversion  = {
            ID=id; PeptideID=peptideid; ModificationID=modificationid; Location=location; Residue=residue; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createModLocationWith id peptideid modificationid location residue rowversion (cvParams:seq<CvParam>) = {
            ID=id; PeptideID=peptideid; ModificationID=modificationid; Location=location; Residue=residue; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (modlocation:ModLocation)  = 
            ParamContainer.addOrUpdateInPlace param modlocation.ParamContainer |> ignore
            modlocation

        let createModLocationTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ModLocation (
                ID INTEGER NOT NULL,
                PeptideID INTEGER NOT NULL,
                ModificationID INTEGER NOT NULL,
                Location INTEGER NOT NULL,
                Residue TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT fk_Peptides_id FOREIGN KEY (PeptideID) REFERENCES Peptide (ID),
                CONSTRAINT fk_Modification_id FOREIGN KEY (ModificationID) REFERENCES Modification (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createModLocationParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ModLocationParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertModLocation (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ModLocation (
                    ID,
                    PeptideID,
                    ModificationID,
                    Location,
                    Residue,
                    RowVersion)
                    VALUES (
                        @id,
                        @peptideid,
                        @modificationid,
                        @location,
                        @residue,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@peptideid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@modificationid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@location",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@residue",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id peptideid modificationid location residue rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@peptideid"].Value <- peptideid
                cmd.Parameters.["@modificationid"].Value <- modificationid
                cmd.Parameters.["@location"].Value <- location
                cmd.Parameters.["@residue"].Value <- residue
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertModLocationParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ModLocationParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertModLocationToDb insertInTableF insertInParamTableF modlocation =
             insertInTableF
                modlocation.ID
                modlocation.PeptideID
                modlocation.ModificationID
                modlocation.Location
                modlocation.Residue
                modlocation.RowVersion
             insertInParamTableF modlocation.ParamContainer

            ///
        let prepareSelectModLocationbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModLocationbyPeptideID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE PeptideID=@peptideid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ peptideid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun peptideid->
            cmd.Parameters.["@ peptideid "].Value <- peptideid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModLocationbyModificationID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE ModificationID=@modificationid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ modificationid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun modificationid->
            cmd.Parameters.["@ modificationid "].Value <- modificationid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModLocationbyLocation (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE Location=@location" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ location ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun location->
            cmd.Parameters.["@ location "].Value <- location
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModLocationbyResidue (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE Residue=@residue" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ residue ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun residue->
            cmd.Parameters.["@ residue "].Value <- residue
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModLocationbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModLocation WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetInt32(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Modification =
        open DataModel
        open ParamContainer


        type Modification = {
            ID : int 
            Name : string option;
            Residues : string option;
            MonoisotopicMassDelta : float option;
            AvgMassDelta : float option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createModification id name residues monoisotopicmassdelta avgmassdelta rowversion  = {
            ID=id; Name=name; Residues=residues; MonoisotopicMassDelta=monoisotopicmassdelta; AvgMassDelta=avgmassdelta; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createModificationWith id name residues monoisotopicmassdelta avgmassdelta rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; Residues=residues; MonoisotopicMassDelta=monoisotopicmassdelta; AvgMassDelta=avgmassdelta; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (modification:Modification)  = 
            ParamContainer.addOrUpdateInPlace param modification.ParamContainer |> ignore
            modification

        let createModificationTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Modification (
                ID INTEGER NOT NULL,
                Name TEXT,
                Residues TEXT,
                MonoisotopicMassDelta REAL,
                AvgMassDelta REAL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createModificationParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ModificationParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertModification (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Modification (
                    ID,
                    Name,
                    Residues,
                    MonoisotopicMassDelta,
                    AvgMassDelta,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @residues,
                        @monoisotopicmassdelta,
                        @avgmassdelta,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@residues",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@monoisotopicmassdelta",Data.DbType.Double) |> ignore
            cmd.Parameters.Add("@avgmassdelta",Data.DbType.Double) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name residues monoisotopicmassdelta avgmassdelta rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@residues"].Value <- residues
                cmd.Parameters.["@monoisotopicmassdelta"].Value <- monoisotopicmassdelta
                cmd.Parameters.["@avgmassdelta"].Value <- avgmassdelta
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertModificationParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ModificationParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertModificationToDb insertInTableF insertInParamTableF modification =
             insertInTableF
                modification.ID
                modification.Name
                modification.Residues
                modification.MonoisotopicMassDelta
                modification.AvgMassDelta
                modification.RowVersion
             insertInParamTableF modification.ParamContainer

            ///
        let prepareSelectModificationbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModificationbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModificationbyResidues (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE Residues=@residues" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ residues ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun residues->
            cmd.Parameters.["@ residues "].Value <- residues
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModificationbyMonoisotopicMassDelta (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE MonoisotopicMassDelta=@monoisotopicmassdelta" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ monoisotopicmassdelta ",Data.DbType.Double) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun monoisotopicmassdelta->
            cmd.Parameters.["@ monoisotopicmassdelta "].Value <- monoisotopicmassdelta
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModificationbyAvgMassDelta (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE AvgMassDelta=@avgmassdelta" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ avgmassdelta ",Data.DbType.Double) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun avgmassdelta->
            cmd.Parameters.["@ avgmassdelta "].Value <- avgmassdelta
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectModificationbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMModification WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDouble(3), reader.GetDouble(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Ontology =
        open DataModel
        open ParamContainer


        type Ontology = {
            ID : int 
            Name : string 
            ParamContainer : ParamContainer
            }

        let createOntology id name  = {
            ID=id; Name=name; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createOntologyWith id name (cvParams:seq<CvParam>) = {
            ID=id; Name=name; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (ontology:Ontology)  = 
            ParamContainer.addOrUpdateInPlace param ontology.ParamContainer |> ignore
            ontology

        let createOntologyTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Ontology (
                ID INTEGER NOT NULL,
                Name TEXT NOT NULL,
                RowVersion  BLOB(8) NOT NULL DEFAULT 0
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createOntologyParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE OntologyParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertOntology (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Ontology (
                    ID,
                    Name)
                    VALUES (
                        @id,
                        @name)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            (fun id name ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.ExecuteNonQuery()
                )

        let prepareInsertOntologyParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO OntologyParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertOntologyToDb insertInTableF insertInParamTableF ontology =
             insertInTableF
                ontology.ID
                ontology.Name
             insertInParamTableF ontology.ParamContainer

            ///
        let prepareSelectOntologybyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOntology WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectOntologybyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOntology WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module TermTag =
        open DataModel
        open ParamContainer


        type TermTag = {
            ID : int 
            TermID : int option;
            Name : string 
            Value : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createTermTag id termid name value rowversion  = {
            ID=id; TermID=termid; Name=name; Value=value; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createTermTagWith id termid name value rowversion (cvParams:seq<CvParam>) = {
            ID=id; TermID=termid; Name=name; Value=value; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (termtag:TermTag)  = 
            ParamContainer.addOrUpdateInPlace param termtag.ParamContainer |> ignore
            termtag

        let createTermTagTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE TermTag (
                ID INTEGER NOT NULL,
                TermID INTEGER,
                Name TEXT NOT NULL,
                Value TEXT NOT NULL,
                RowVersion  BLOB(8) NOT NULL DEFAULT 0,
                CONSTRAINT FK_TermTag_Term FOREIGN KEY (TermID) REFERENCES Term (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createTermTagParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE TermTagParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertTermTag (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO TermTag (
                    ID,
                    TermID,
                    Name,
                    Value,
                    RowVersion)
                    VALUES (
                        @id,
                        @termid,
                        @name,
                        @value,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@termid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id termid name value rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@termid"].Value <- termid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertTermTagParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO TermTagParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertTermTagToDb insertInTableF insertInParamTableF termtag =
             insertInTableF
                termtag.ID
                termtag.TermID
                termtag.Name
                termtag.Value
                termtag.RowVersion
             insertInParamTableF termtag.ParamContainer

            ///
        let prepareSelectTermTagbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermTag WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermTagbyTermID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermTag WHERE TermID=@termid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ termid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun termid->
            cmd.Parameters.["@ termid "].Value <- termid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermTagbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermTag WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermTagbyValue (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermTag WHERE Value=@value" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ value ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun value->
            cmd.Parameters.["@ value "].Value <- value
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermTagbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermTag WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module TermRelationShip =
        open DataModel
        open ParamContainer


        type TermRelationShip = {
            ID : int 
            TermID : int option;
            RelationShipType : string 
            FKRelatedTerm : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createTermRelationShip id termid relationshiptype fkrelatedterm rowversion  = {
            ID=id; TermID=termid; RelationShipType=relationshiptype; FKRelatedTerm=fkrelatedterm; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createTermRelationShipWith id termid relationshiptype fkrelatedterm rowversion (cvParams:seq<CvParam>) = {
            ID=id; TermID=termid; RelationShipType=relationshiptype; FKRelatedTerm=fkrelatedterm; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (termrelationship:TermRelationShip)  = 
            ParamContainer.addOrUpdateInPlace param termrelationship.ParamContainer |> ignore
            termrelationship

        let createTermRelationShipTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE TermRelationShip (
                ID INTEGER NOT NULL,
                TermID INTEGER,
                RelationShipType TEXT NOT NULL,
                FKRelatedTerm TEXT NOT NULL,
                RowVersion  BLOB(8) NOT NULL DEFAULT 0,
                CONSTRAINT FK_TermRelationShip_Term_ID FOREIGN KEY (TermID) REFERENCES Term (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createTermRelationShipParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE TermRelationShipParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertTermRelationShip (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO TermRelationShip (
                    ID,
                    TermID,
                    RelationShipType,
                    FKRelatedTerm,
                    RowVersion)
                    VALUES (
                        @id,
                        @termid,
                        @relationshiptype,
                        @fkrelatedterm,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@termid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@relationshiptype",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@fkrelatedterm",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id termid relationshiptype fkrelatedterm rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@termid"].Value <- termid
                cmd.Parameters.["@relationshiptype"].Value <- relationshiptype
                cmd.Parameters.["@fkrelatedterm"].Value <- fkrelatedterm
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertTermRelationShipParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO TermRelationShipParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertTermRelationShipToDb insertInTableF insertInParamTableF termrelationship =
             insertInTableF
                termrelationship.ID
                termrelationship.TermID
                termrelationship.RelationShipType
                termrelationship.FKRelatedTerm
                termrelationship.RowVersion
             insertInParamTableF termrelationship.ParamContainer

            ///
        let prepareSelectTermRelationShipbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermRelationShip WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermRelationShipbyTermID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermRelationShip WHERE TermID=@termid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ termid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun termid->
            cmd.Parameters.["@ termid "].Value <- termid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermRelationShipbyRelationShipType (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermRelationShip WHERE RelationShipType=@relationshiptype" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ relationshiptype ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun relationshiptype->
            cmd.Parameters.["@ relationshiptype "].Value <- relationshiptype
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermRelationShipbyFKRelatedTerm (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermRelationShip WHERE FKRelatedTerm=@fkrelatedterm" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ fkrelatedterm ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun fkrelatedterm->
            cmd.Parameters.["@ fkrelatedterm "].Value <- fkrelatedterm
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectTermRelationShipbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMTermRelationShip WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ProteinDetectionHypothesis =
        open DataModel
        open ParamContainer


        type ProteinDetectionHypothesis = {
            ID : int 
            DBSequenceID : int 
            ProteinAmbiguityGroupID : int 
            Name : string option;
            PassThreshold : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createProteinDetectionHypothesis id dbsequenceid proteinambiguitygroupid name passthreshold rowversion  = {
            ID=id; DBSequenceID=dbsequenceid; ProteinAmbiguityGroupID=proteinambiguitygroupid; Name=name; PassThreshold=passthreshold; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createProteinDetectionHypothesisWith id dbsequenceid proteinambiguitygroupid name passthreshold rowversion (cvParams:seq<CvParam>) = {
            ID=id; DBSequenceID=dbsequenceid; ProteinAmbiguityGroupID=proteinambiguitygroupid; Name=name; PassThreshold=passthreshold; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (proteindetectionhypothesis:ProteinDetectionHypothesis)  = 
            ParamContainer.addOrUpdateInPlace param proteindetectionhypothesis.ParamContainer |> ignore
            proteindetectionhypothesis

        let createProteinDetectionHypothesisTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionHypothesis (
                ID INTEGER NOT NULL,
                DBSequenceID INTEGER NOT NULL,
                ProteinAmbiguityGroupID INTEGER NOT NULL,
                Name TEXT,
                PassThreshold TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_ProtDetectHypothesis_ProtAmbigGroup_ID FOREIGN KEY (ProteinAmbiguityGroupID) REFERENCES ProteinAmbiguityGroup (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createProteinDetectionHypothesisParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionHypothesisParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertProteinDetectionHypothesis (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionHypothesis (
                    ID,
                    DBSequenceID,
                    ProteinAmbiguityGroupID,
                    Name,
                    PassThreshold,
                    RowVersion)
                    VALUES (
                        @id,
                        @dbsequenceid,
                        @proteinambiguitygroupid,
                        @name,
                        @passthreshold,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@dbsequenceid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@proteinambiguitygroupid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@passthreshold",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id dbsequenceid proteinambiguitygroupid name passthreshold rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@dbsequenceid"].Value <- dbsequenceid
                cmd.Parameters.["@proteinambiguitygroupid"].Value <- proteinambiguitygroupid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@passthreshold"].Value <- passthreshold
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertProteinDetectionHypothesisParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionHypothesisParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertProteinDetectionHypothesisToDb insertInTableF insertInParamTableF proteindetectionhypothesis =
             insertInTableF
                proteindetectionhypothesis.ID
                proteindetectionhypothesis.DBSequenceID
                proteindetectionhypothesis.ProteinAmbiguityGroupID
                proteindetectionhypothesis.Name
                proteindetectionhypothesis.PassThreshold
                proteindetectionhypothesis.RowVersion
             insertInParamTableF proteindetectionhypothesis.ParamContainer

            ///
        let prepareSelectProteinDetectionHypothesisbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionHypothesisbyDBSequenceID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE DBSequenceID=@dbsequenceid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ dbsequenceid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun dbsequenceid->
            cmd.Parameters.["@ dbsequenceid "].Value <- dbsequenceid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionHypothesisbyProteinAmbiguityGroupID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE ProteinAmbiguityGroupID=@proteinambiguitygroupid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ proteinambiguitygroupid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun proteinambiguitygroupid->
            cmd.Parameters.["@ proteinambiguitygroupid "].Value <- proteinambiguitygroupid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionHypothesisbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionHypothesisbyPassThreshold (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE PassThreshold=@passthreshold" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ passthreshold ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun passthreshold->
            cmd.Parameters.["@ passthreshold "].Value <- passthreshold
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionHypothesisbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionHypothesis WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetInt32(2), reader.GetString(3), reader.GetString(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ProteinDetectionList =
        open DataModel
        open ParamContainer


        type ProteinDetectionList = {
            ID : int 
            Accession : string 
            Name : string 
            SearchDBID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createProteinDetectionList id accession name searchdbid rowversion  = {
            ID=id; Accession=accession; Name=name; SearchDBID=searchdbid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createProteinDetectionListWith id accession name searchdbid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Accession=accession; Name=name; SearchDBID=searchdbid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (proteindetectionlist:ProteinDetectionList)  = 
            ParamContainer.addOrUpdateInPlace param proteindetectionlist.ParamContainer |> ignore
            proteindetectionlist

        let createProteinDetectionListTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionList (
                ID INTEGER NOT NULL,
                Accession TEXT NOT NULL,
                Name TEXT NOT NULL,
                SearchDBID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createProteinDetectionListParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionListParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertProteinDetectionList (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionList (
                    ID,
                    Accession,
                    Name,
                    SearchDBID,
                    RowVersion)
                    VALUES (
                        @id,
                        @accession,
                        @name,
                        @searchdbid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@accession",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@searchdbid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id accession name searchdbid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@accession"].Value <- accession
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@searchdbid"].Value <- searchdbid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertProteinDetectionListParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionListParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertProteinDetectionListToDb insertInTableF insertInParamTableF proteindetectionlist =
             insertInTableF
                proteindetectionlist.ID
                proteindetectionlist.Accession
                proteindetectionlist.Name
                proteindetectionlist.SearchDBID
                proteindetectionlist.RowVersion
             insertInParamTableF proteindetectionlist.ParamContainer

            ///
        let prepareSelectProteinDetectionListbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionList WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionListbyAccession (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionList WHERE Accession=@accession" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ accession ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun accession->
            cmd.Parameters.["@ accession "].Value <- accession
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionListbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionList WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionListbySearchDBID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionList WHERE SearchDBID=@searchdbid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ searchdbid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun searchdbid->
            cmd.Parameters.["@ searchdbid "].Value <- searchdbid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionListbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionList WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetString(3), reader.GetDateTime(4) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module PeptideHypothesis =
        open DataModel
        open ParamContainer


        type PeptideHypothesis = {
            ID : int 
            PeptideEvidenceID : string 
            ProteinDetectionHypothesisID : int option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createPeptideHypothesis id peptideevidenceid proteindetectionhypothesisid rowversion  = {
            ID=id; PeptideEvidenceID=peptideevidenceid; ProteinDetectionHypothesisID=proteindetectionhypothesisid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createPeptideHypothesisWith id peptideevidenceid proteindetectionhypothesisid rowversion (cvParams:seq<CvParam>) = {
            ID=id; PeptideEvidenceID=peptideevidenceid; ProteinDetectionHypothesisID=proteindetectionhypothesisid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (peptidehypothesis:PeptideHypothesis)  = 
            ParamContainer.addOrUpdateInPlace param peptidehypothesis.ParamContainer |> ignore
            peptidehypothesis

        let createPeptideHypothesisTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PeptideHypothesis (
                ID INTEGER NOT NULL,
                PeptideEvidenceID TEXT NOT NULL,
                ProteinDetectionHypothesisID INTEGER,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_PepHypothesis_PepEvidence_ID FOREIGN KEY (PeptideEvidenceID) REFERENCES PeptideEvidence (ID),
                CONSTRAINT FK_PepHypothesis_ProtDetectHypo_ID FOREIGN KEY (ProteinDetectionHypothesisID) REFERENCES ProteinDetectionHypothesis (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createPeptideHypothesisParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PeptideHypothesisParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertPeptideHypothesis (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PeptideHypothesis (
                    ID,
                    PeptideEvidenceID,
                    ProteinDetectionHypothesisID,
                    RowVersion)
                    VALUES (
                        @id,
                        @peptideevidenceid,
                        @proteindetectionhypothesisid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@peptideevidenceid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@proteindetectionhypothesisid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id peptideevidenceid proteindetectionhypothesisid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@peptideevidenceid"].Value <- peptideevidenceid
                cmd.Parameters.["@proteindetectionhypothesisid"].Value <- proteindetectionhypothesisid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertPeptideHypothesisParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PeptideHypothesisParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertPeptideHypothesisToDb insertInTableF insertInParamTableF peptidehypothesis =
             insertInTableF
                peptidehypothesis.ID
                peptidehypothesis.PeptideEvidenceID
                peptidehypothesis.ProteinDetectionHypothesisID
                peptidehypothesis.RowVersion
             insertInParamTableF peptidehypothesis.ParamContainer

            ///
        let prepareSelectPeptideHypothesisbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideHypothesis WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideHypothesisbyPeptideEvidenceID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideHypothesis WHERE PeptideEvidenceID=@peptideevidenceid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ peptideevidenceid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun peptideevidenceid->
            cmd.Parameters.["@ peptideevidenceid "].Value <- peptideevidenceid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideHypothesisbyProteinDetectionHypothesisID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideHypothesis WHERE ProteinDetectionHypothesisID=@proteindetectionhypothesisid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ proteindetectionhypothesisid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun proteindetectionhypothesisid->
            cmd.Parameters.["@ proteindetectionhypothesisid "].Value <- proteindetectionhypothesisid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPeptideHypothesisbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPeptideHypothesis WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetInt32(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ProteinAmbiguityGroup =
        open DataModel
        open ParamContainer


        type ProteinAmbiguityGroup = {
            ID : int 
            ProteinDetectionListID : int 
            Name : string option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createProteinAmbiguityGroup id proteindetectionlistid name rowversion  = {
            ID=id; ProteinDetectionListID=proteindetectionlistid; Name=name; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createProteinAmbiguityGroupWith id proteindetectionlistid name rowversion (cvParams:seq<CvParam>) = {
            ID=id; ProteinDetectionListID=proteindetectionlistid; Name=name; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (proteinambiguitygroup:ProteinAmbiguityGroup)  = 
            ParamContainer.addOrUpdateInPlace param proteinambiguitygroup.ParamContainer |> ignore
            proteinambiguitygroup

        let createProteinAmbiguityGroupTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinAmbiguityGroup (
                ID INTEGER NOT NULL,
                ProteinDetectionListID INTEGER NOT NULL,
                Name TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_ProtAmbigGroup_ProtDetectList_ID FOREIGN KEY (ProteinDetectionListID) REFERENCES ProteinDetectionList (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createProteinAmbiguityGroupParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinAmbiguityGroupParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertProteinAmbiguityGroup (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinAmbiguityGroup (
                    ID,
                    ProteinDetectionListID,
                    Name,
                    RowVersion)
                    VALUES (
                        @id,
                        @proteindetectionlistid,
                        @name,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@proteindetectionlistid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id proteindetectionlistid name rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@proteindetectionlistid"].Value <- proteindetectionlistid
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertProteinAmbiguityGroupParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinAmbiguityGroupParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertProteinAmbiguityGroupToDb insertInTableF insertInParamTableF proteinambiguitygroup =
             insertInTableF
                proteinambiguitygroup.ID
                proteinambiguitygroup.ProteinDetectionListID
                proteinambiguitygroup.Name
                proteinambiguitygroup.RowVersion
             insertInParamTableF proteinambiguitygroup.ParamContainer

            ///
        let prepareSelectProteinAmbiguityGroupbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinAmbiguityGroup WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinAmbiguityGroupbyProteinDetectionListID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinAmbiguityGroup WHERE ProteinDetectionListID=@proteindetectionlistid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ proteindetectionlistid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun proteindetectionlistid->
            cmd.Parameters.["@ proteindetectionlistid "].Value <- proteindetectionlistid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinAmbiguityGroupbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinAmbiguityGroup WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinAmbiguityGroupbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinAmbiguityGroup WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetInt32(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ProteinDetection =
        open DataModel
        open ParamContainer


        type ProteinDetection = {
            ID : int 
            Name : string option;
            ActivityDate : string option;
            ProteinDetectionListID : int 
            ProteinDetectionProtocolID : int 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createProteinDetection id name activitydate proteindetectionlistid proteindetectionprotocolid rowversion  = {
            ID=id; Name=name; ActivityDate=activitydate; ProteinDetectionListID=proteindetectionlistid; ProteinDetectionProtocolID=proteindetectionprotocolid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createProteinDetectionWith id name activitydate proteindetectionlistid proteindetectionprotocolid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; ActivityDate=activitydate; ProteinDetectionListID=proteindetectionlistid; ProteinDetectionProtocolID=proteindetectionprotocolid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (proteindetection:ProteinDetection)  = 
            ParamContainer.addOrUpdateInPlace param proteindetection.ParamContainer |> ignore
            proteindetection

        let createProteinDetectionTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetection (
                ID INTEGER NOT NULL,
                Name TEXT,
                ActivityDate TEXT,
                ProteinDetectionListID INTEGER NOT NULL,
                ProteinDetectionProtocolID INTEGER NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_ProtDetection_ProtDetecionList_ID FOREIGN KEY (ProteinDetectionListID) REFERENCES ProteinDetectionList (ID),
                CONSTRAINT FK_ProtDetection_ProteinDetectionProtocol_ID FOREIGN KEY (ProteinDetectionProtocolID) REFERENCES ProteinDetectionProtocol (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createProteinDetectionParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertProteinDetection (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetection (
                    ID,
                    Name,
                    ActivityDate,
                    ProteinDetectionListID,
                    ProteinDetectionProtocolID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @activitydate,
                        @proteindetectionlistid,
                        @proteindetectionprotocolid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@activitydate",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@proteindetectionlistid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@proteindetectionprotocolid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name activitydate proteindetectionlistid proteindetectionprotocolid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@activitydate"].Value <- activitydate
                cmd.Parameters.["@proteindetectionlistid"].Value <- proteindetectionlistid
                cmd.Parameters.["@proteindetectionprotocolid"].Value <- proteindetectionprotocolid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertProteinDetectionParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertProteinDetectionToDb insertInTableF insertInParamTableF proteindetection =
             insertInTableF
                proteindetection.ID
                proteindetection.Name
                proteindetection.ActivityDate
                proteindetection.ProteinDetectionListID
                proteindetection.ProteinDetectionProtocolID
                proteindetection.RowVersion
             insertInParamTableF proteindetection.ParamContainer

            ///
        let prepareSelectProteinDetectionbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionbyActivityDate (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE ActivityDate=@activitydate" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ activitydate ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun activitydate->
            cmd.Parameters.["@ activitydate "].Value <- activitydate
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionbyProteinDetectionListID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE ProteinDetectionListID=@proteindetectionlistid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ proteindetectionlistid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun proteindetectionlistid->
            cmd.Parameters.["@ proteindetectionlistid "].Value <- proteindetectionlistid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionbyProteinDetectionProtocolID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE ProteinDetectionProtocolID=@proteindetectionprotocolid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ proteindetectionprotocolid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun proteindetectionprotocolid->
            cmd.Parameters.["@ proteindetectionprotocolid "].Value <- proteindetectionprotocolid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetection WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module SpectrumIdentification =
        open DataModel
        open ParamContainer


        type SpectrumIdentification = {
            ID : int 
            Name : string option;
            ActivityDate : string option;
            SpectrumIdentficationListID : int 
            SpectrumIdentficationProtocolID : int 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createSpectrumIdentification id name activitydate spectrumidentficationlistid spectrumidentficationprotocolid rowversion  = {
            ID=id; Name=name; ActivityDate=activitydate; SpectrumIdentficationListID=spectrumidentficationlistid; SpectrumIdentficationProtocolID=spectrumidentficationprotocolid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createSpectrumIdentificationWith id name activitydate spectrumidentficationlistid spectrumidentficationprotocolid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; ActivityDate=activitydate; SpectrumIdentficationListID=spectrumidentficationlistid; SpectrumIdentficationProtocolID=spectrumidentficationprotocolid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (spectrumidentification:SpectrumIdentification)  = 
            ParamContainer.addOrUpdateInPlace param spectrumidentification.ParamContainer |> ignore
            spectrumidentification

        let createSpectrumIdentificationTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentification (
                ID INTEGER NOT NULL,
                Name TEXT,
                ActivityDate TEXT,
                SpectrumIdentficationListID INTEGER NOT NULL,
                SpectrumIdentficationProtocolID INTEGER NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 8,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_SpectrumIdentification_SpectrumIdentifcationList_ID FOREIGN KEY (SpectrumIdentficationListID) REFERENCES SpectrumIdentificationList (ID),
                CONSTRAINT FK_SpectrumIdentification_SpectrumIdentificationProtocol_ID FOREIGN KEY () REFERENCES SpectrumIdentificationProtocol (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createSpectrumIdentificationParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertSpectrumIdentification (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentification (
                    ID,
                    Name,
                    ActivityDate,
                    SpectrumIdentficationListID,
                    SpectrumIdentficationProtocolID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @activitydate,
                        @spectrumidentficationlistid,
                        @spectrumidentficationprotocolid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@activitydate",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@spectrumidentficationlistid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@spectrumidentficationprotocolid",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name activitydate spectrumidentficationlistid spectrumidentficationprotocolid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@activitydate"].Value <- activitydate
                cmd.Parameters.["@spectrumidentficationlistid"].Value <- spectrumidentficationlistid
                cmd.Parameters.["@spectrumidentficationprotocolid"].Value <- spectrumidentficationprotocolid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertSpectrumIdentificationParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertSpectrumIdentificationToDb insertInTableF insertInParamTableF spectrumidentification =
             insertInTableF
                spectrumidentification.ID
                spectrumidentification.Name
                spectrumidentification.ActivityDate
                spectrumidentification.SpectrumIdentficationListID
                spectrumidentification.SpectrumIdentficationProtocolID
                spectrumidentification.RowVersion
             insertInParamTableF spectrumidentification.ParamContainer

            ///
        let prepareSelectSpectrumIdentificationbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationbyActivityDate (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE ActivityDate=@activitydate" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ activitydate ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun activitydate->
            cmd.Parameters.["@ activitydate "].Value <- activitydate
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationbySpectrumIdentficationListID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE SpectrumIdentficationListID=@spectrumidentficationlistid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectrumidentficationlistid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun spectrumidentficationlistid->
            cmd.Parameters.["@ spectrumidentficationlistid "].Value <- spectrumidentficationlistid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationbySpectrumIdentficationProtocolID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE SpectrumIdentficationProtocolID=@spectrumidentficationprotocolid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ spectrumidentficationprotocolid ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun spectrumidentficationprotocolid->
            cmd.Parameters.["@ spectrumidentficationprotocolid "].Value <- spectrumidentficationprotocolid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentification WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetInt32(3), reader.GetInt32(4), reader.GetDateTime(5) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module SpectrumIdentificationProtocol =
        open DataModel
        open ParamContainer


        type SpectrumIdentificationProtocol = {
            ID : int 
            Name : string option;
            AnalysisSoftwareID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createSpectrumIdentificationProtocol id name analysissoftwareid rowversion  = {
            ID=id; Name=name; AnalysisSoftwareID=analysissoftwareid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createSpectrumIdentificationProtocolWith id name analysissoftwareid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; AnalysisSoftwareID=analysissoftwareid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (spectrumidentificationprotocol:SpectrumIdentificationProtocol)  = 
            ParamContainer.addOrUpdateInPlace param spectrumidentificationprotocol.ParamContainer |> ignore
            spectrumidentificationprotocol

        let createSpectrumIdentificationProtocolTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationProtocol (
                ID INTEGER NOT NULL,
                Name TEXT,
                AnalysisSoftwareID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK-SpecIdentProtocol_AnalysisSoftware_ID FOREIGN KEY (AnalysisSoftwareID) REFERENCES AnalysisSoftware (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createSpectrumIdentificationProtocolParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE SpectrumIdentificationProtocolParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertSpectrumIdentificationProtocol (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationProtocol (
                    ID,
                    Name,
                    AnalysisSoftwareID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @analysissoftwareid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@analysissoftwareid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name analysissoftwareid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@analysissoftwareid"].Value <- analysissoftwareid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertSpectrumIdentificationProtocolParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO SpectrumIdentificationProtocolParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertSpectrumIdentificationProtocolToDb insertInTableF insertInParamTableF spectrumidentificationprotocol =
             insertInTableF
                spectrumidentificationprotocol.ID
                spectrumidentificationprotocol.Name
                spectrumidentificationprotocol.AnalysisSoftwareID
                spectrumidentificationprotocol.RowVersion
             insertInParamTableF spectrumidentificationprotocol.ParamContainer

            ///
        let prepareSelectSpectrumIdentificationProtocolbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationProtocol WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationProtocolbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationProtocol WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationProtocolbyAnalysisSoftwareID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationProtocol WHERE AnalysisSoftwareID=@analysissoftwareid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ analysissoftwareid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun analysissoftwareid->
            cmd.Parameters.["@ analysissoftwareid "].Value <- analysissoftwareid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectSpectrumIdentificationProtocolbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMSpectrumIdentificationProtocol WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module ProteinDetectionProtocol =
        open DataModel
        open ParamContainer


        type ProteinDetectionProtocol = {
            ID : int 
            Name : string option;
            AnalysisSoftwareID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createProteinDetectionProtocol id name analysissoftwareid rowversion  = {
            ID=id; Name=name; AnalysisSoftwareID=analysissoftwareid; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createProteinDetectionProtocolWith id name analysissoftwareid rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; AnalysisSoftwareID=analysissoftwareid; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (proteindetectionprotocol:ProteinDetectionProtocol)  = 
            ParamContainer.addOrUpdateInPlace param proteindetectionprotocol.ParamContainer |> ignore
            proteindetectionprotocol

        let createProteinDetectionProtocolTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionProtocol (
                ID INTEGER NOT NULL,
                Name TEXT,
                AnalysisSoftwareID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 8,
                PRIMARY KEY (ID) ,
                CONSTRAINT ProteinDetectionProt_AnalysisSoftware_ID FOREIGN KEY (AnalysisSoftwareID) REFERENCES AnalysisSoftware (ID)
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createProteinDetectionProtocolParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE ProteinDetectionProtocolParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertProteinDetectionProtocol (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionProtocol (
                    ID,
                    Name,
                    AnalysisSoftwareID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @analysissoftwareid,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@analysissoftwareid",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name analysissoftwareid rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@analysissoftwareid"].Value <- analysissoftwareid
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertProteinDetectionProtocolParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO ProteinDetectionProtocolParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertProteinDetectionProtocolToDb insertInTableF insertInParamTableF proteindetectionprotocol =
             insertInTableF
                proteindetectionprotocol.ID
                proteindetectionprotocol.Name
                proteindetectionprotocol.AnalysisSoftwareID
                proteindetectionprotocol.RowVersion
             insertInParamTableF proteindetectionprotocol.ParamContainer

            ///
        let prepareSelectProteinDetectionProtocolbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionProtocol WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionProtocolbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionProtocol WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionProtocolbyAnalysisSoftwareID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionProtocol WHERE AnalysisSoftwareID=@analysissoftwareid" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ analysissoftwareid ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun analysissoftwareid->
            cmd.Parameters.["@ analysissoftwareid "].Value <- analysissoftwareid
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectProteinDetectionProtocolbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMProteinDetectionProtocol WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Person =
        open DataModel
        open ParamContainer


        type Person = {
            ID : int 
            Name : string option;
            AnalysisSoftware_ID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createPerson id name analysissoftware_id rowversion  = {
            ID=id; Name=name; AnalysisSoftware_ID=analysissoftware_id; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createPersonWith id name analysissoftware_id rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; AnalysisSoftware_ID=analysissoftware_id; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (person:Person)  = 
            ParamContainer.addOrUpdateInPlace param person.ParamContainer |> ignore
            person

        let createPersonTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Person (
                ID INTEGER NOT NULL,
                Name TEXT,
                AnalysisSoftware_ID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createPersonParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE PersonParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertPerson (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Person (
                    ID,
                    Name,
                    AnalysisSoftware_ID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @analysissoftware_id,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@analysissoftware_id",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name analysissoftware_id rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@analysissoftware_id"].Value <- analysissoftware_id
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertPersonParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO PersonParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertPersonToDb insertInTableF insertInParamTableF person =
             insertInTableF
                person.ID
                person.Name
                person.AnalysisSoftware_ID
                person.RowVersion
             insertInParamTableF person.ParamContainer

            ///
        let prepareSelectPersonbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPerson WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPersonbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPerson WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPersonbyAnalysisSoftware_ID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPerson WHERE AnalysisSoftware_ID=@analysissoftware_id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ analysissoftware_id ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun analysissoftware_id->
            cmd.Parameters.["@ analysissoftware_id "].Value <- analysissoftware_id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectPersonbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMPerson WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module Organization =
        open DataModel
        open ParamContainer


        type Organization = {
            ID : int 
            Name : string option;
            AnalysisSoftware_ID : string 
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createOrganization id name analysissoftware_id rowversion  = {
            ID=id; Name=name; AnalysisSoftware_ID=analysissoftware_id; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createOrganizationWith id name analysissoftware_id rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; AnalysisSoftware_ID=analysissoftware_id; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (organization:Organization)  = 
            ParamContainer.addOrUpdateInPlace param organization.ParamContainer |> ignore
            organization

        let createOrganizationTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE Organization (
                ID INTEGER NOT NULL,
                Name TEXT,
                AnalysisSoftware_ID TEXT NOT NULL,
                RowVersion BLOB(8) NOT NULL DEFAULT 8,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createOrganizationParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE OrganizationParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertOrganization (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO Organization (
                    ID,
                    Name,
                    AnalysisSoftware_ID,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @analysissoftware_id,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@analysissoftware_id",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name analysissoftware_id rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@analysissoftware_id"].Value <- analysissoftware_id
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertOrganizationParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO OrganizationParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertOrganizationToDb insertInTableF insertInParamTableF organization =
             insertInTableF
                organization.ID
                organization.Name
                organization.AnalysisSoftware_ID
                organization.RowVersion
             insertInParamTableF organization.ParamContainer

            ///
        let prepareSelectOrganizationbyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOrganization WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectOrganizationbyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOrganization WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectOrganizationbyAnalysisSoftware_ID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOrganization WHERE AnalysisSoftware_ID=@analysissoftware_id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ analysissoftware_id ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun analysissoftware_id->
            cmd.Parameters.["@ analysissoftware_id "].Value <- analysissoftware_id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectOrganizationbyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMOrganization WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetString(2), reader.GetDateTime(3) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []

    module AnalysisSoftware =
        open DataModel
        open ParamContainer


        type AnalysisSoftware = {
            ID : int 
            Name : string option;
            RowVersion : System.DateTime 
            ParamContainer : ParamContainer
            }

        let createAnalysisSoftware id name rowversion  = {
            ID=id; Name=name; RowVersion=rowversion; ParamContainer=System.Collections.Generic.Dictionary<TermId,CvParam>()}

        let createAnalysisSoftwareWith id name rowversion (cvParams:seq<CvParam>) = {
            ID=id; Name=name; RowVersion=rowversion; ParamContainer=ParamContainer.ofSeq cvParams}

        let addOrUpDateInPlace (param:CvParam)  (analysissoftware:AnalysisSoftware)  = 
            ParamContainer.addOrUpdateInPlace param analysissoftware.ParamContainer |> ignore
            analysissoftware

        let createAnalysisSoftwareTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE AnalysisSoftware (
                ID INTEGER NOT NULL,
                Name TEXT,
                RowVersion BLOB(8) NOT NULL,
                PRIMARY KEY (ID) 
                )" 
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let createAnalysisSoftwareParamTable (cn:SQLiteConnection) =
            let querystring = 
                "CREATE TABLE AnalysisSoftwareParam (
                ID INTEGER NOT NULL,
                FKParamContainer INTEGER NOT NULL,
                FKTerm INTEGER NOT NULL,
                FKUnit INTEGER,
                Value TEXT,
                RowVersion BLOB(8) NOT NULL DEFAULT 0,
                PRIMARY KEY (ID) ,
                CONSTRAINT FK_Param_ProteinDetectionProtocol_ID FOREIGN KEY (FKParamContainer) REFERENCES ProteinDetectionProtocol (ID),
                CONSTRAINT FK_ProtDetectProtocolTerm_Term_ID FOREIGN KEY (FKTerm) REFERENCES Term (ID),
                CONSTRAINT FK_ProtDetectProtocolUnit_Term_ID FOREIGN KEY (FKUnit) REFERENCES Term (ID))"
            let cmd  = new SQLiteCommand(querystring, cn)
            cmd.ExecuteNonQuery()

        let prepareInsertAnalysisSoftware (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO AnalysisSoftware (
                    ID,
                    Name,
                    RowVersion)
                    VALUES (
                        @id,
                        @name,
                        @rowversion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@name",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowversion",Data.DbType.DateTime) |> ignore
            (fun id name rowversion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@name"].Value <- name
                cmd.Parameters.["@rowversion"].Value <- rowversion
                cmd.ExecuteNonQuery()
                )

        let prepareInsertAnalysisSoftwareParam (cn:SQLiteConnection) tr =
            let querystring = 
                "INSERT INTO AnalysisSoftwareParam (
                    ID,
                    FKParamContainer,
                    FKTerm,
                    FKUnit,
                    Value,
                    RowVersion,
                    VALUES (
                        @id,
                        @fkParamContainer,
                        @fkTerm,
                        @fkUnit,
                        @value,
                        @rowVersion)"
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@id",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkParamContainer",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkTerm",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@fkUnit",Data.DbType.Int32) |> ignore
            cmd.Parameters.Add("@value",Data.DbType.String) |> ignore
            cmd.Parameters.Add("@rowVersion",Data.DbType.DateTime) |> ignore
            (fun id fkParamContainer fkTerm fkUnit value rowVersion ->
                cmd.Parameters.["@id"].Value <- id
                cmd.Parameters.["@fkParamContainer"].Value <- fkParamContainer
                cmd.Parameters.["@fkTerm"].Value <- fkTerm
                cmd.Parameters.["@fkUnit"].Value <- fkUnit
                cmd.Parameters.["@value"].Value <- value
                cmd.Parameters.["@rowVersion"].Value <- rowVersion
                cmd.ExecuteNonQuery()
                )

        let insertAnalysisSoftwareToDb insertInTableF insertInParamTableF analysissoftware =
             insertInTableF
                analysissoftware.ID
                analysissoftware.Name
                analysissoftware.RowVersion
             insertInParamTableF analysissoftware.ParamContainer

            ///
        let prepareSelectAnalysisSoftwarebyID (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMAnalysisSoftware WHERE ID=@id" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ id ",Data.DbType.Int32) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun id->
            cmd.Parameters.["@ id "].Value <- id
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectAnalysisSoftwarebyName (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMAnalysisSoftware WHERE Name=@name" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ name ",Data.DbType.String) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun name->
            cmd.Parameters.["@ name "].Value <- name
            use reader = cmd.ExecuteReader() 
            readerloop reader []

            ///
        let prepareSelectAnalysisSoftwarebyRowVersion (cn:SQLiteConnection) tr =
            let querystring = 
                "SELECT * FROMAnalysisSoftware WHERE RowVersion=@rowversion" 
            let cmd = new SQLiteCommand(querystring, cn, tr)
            cmd.Parameters.Add("@ rowversion ",Data.DbType.DateTime) |> ignore
            let rec readerloop (reader:SQLiteDataReader) (acc) =
                match reader.Read() with
                | true  -> readerloop reader  ((reader.GetInt32(0), reader.GetString(1), reader.GetDateTime(2) ) :: acc)
                | false ->  acc 
            fun rowversion->
            cmd.Parameters.["@ rowversion "].Value <- rowversion
            use reader = cmd.ExecuteReader() 
            readerloop reader []


    open FSharp.Care.Monads.Either

    let initDB fileName =
        let _ = FSharp.Care.IO.FileIO.DeleteFile fileName 
        let connectionString = sprintf "Data Source=%s;Version=3" fileName
        use cn = new SQLiteConnection(connectionString)
        cn.Open()   
        DBSequence.createDBSequenceTable      |> ignore
        DBSequence.createDBSequenceParamTable |> ignore
  