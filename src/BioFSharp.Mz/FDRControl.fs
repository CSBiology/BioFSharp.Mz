namespace BioFSharp.Mz

open FSharp.Stats

module FDRControl = 

    /// for given data, creates a logistic regression model and returns a mapping function for this model
    let private getLogisticRegressionFunction (x:vector) (y:vector) epsilon = 
        let alpha = 
            match FSharp.Stats.Fitting.LogisticRegression.Univariable.estimateAlpha epsilon x y with 
            | Some a -> a
            | None -> failwith "Could not find an alpha for logistic regression of fdr data"
        let weight = FSharp.Stats.Fitting.LogisticRegression.Univariable.coefficient epsilon alpha x y
        FSharp.Stats.Fitting.LogisticRegression.Univariable.fit weight

    /// returns scores, pep, q
    let private binningFunction bandwidth pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[])  =   
        data
        |> Array.groupBy (fun s -> floor (scoreF s / bandwidth))
        |> Array.sortBy fst
        |> Array.map (fun (k,values)->
                        let median     = values |> Array.map scoreF |> Array.average
                        let totalCount = values |> Array.length |> float
                        let decoyCount = values |> Array.filter isDecoyF |> Array.length |> float
                        //(median |> float,(decoyCount * pi0  / totalCount))
                        median,totalCount,decoyCount
                            //(median, totalCount )
                     )
        |> fun a ->
            a
            |> Array.mapi (fun i (median,totalCountBin,decoyCountBin) ->
                            /// TODO: Accumulate totalCount + totalDecoyCount beforeHand and skip the time intensive mapping accross the array in each iteration.
                            let _,totalCountRight,decoyCountRight = a.[i..a.Length-1] |> Array.reduce (fun (x,y,z) (x',y',z') -> x+x',y+y',z+z')
                            (median,(pi0 * 2. * decoyCountBin / totalCountBin),(pi0 * 2. * decoyCountRight / totalCountRight))
                          )
        |> Array.sortBy (fun (score,pep,q) -> score) 
        |> Array.unzip3 
        |> fun (score,pep,q) -> vector score, vector pep, vector q
    
    /// Calculates q value mapping funtion for target/decoy dataset
    let getQValueFunc pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let bw = 
            data 
            |> Array.map scoreF
            |> FSharp.Stats.Distributions.Bandwidth.nrd0
        let (scores,_,q) = binningFunction bw pi0 scoreF isDecoyF data
        getLogisticRegressionFunction scores q 0.0000001

    /// Calculates q values for target/decoy dataset
    let getQValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let f = getQValueFunc pi0 scoreF isDecoyF data
        Array.map (scoreF >> f) data

    /// Calculates pep value mapping funtion for target/decoy dataset
    let getPEPValueFunc pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let bw = 
            data 
            |> Array.map scoreF
            |> FSharp.Stats.Distributions.Bandwidth.nrd0
        let (scores,pep,_) = binningFunction bw pi0 scoreF isDecoyF data
        getLogisticRegressionFunction scores pep 0.0000001

    /// Calculates pep values for target/decoy dataset
    let getPEPValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) = 
        let f = getPEPValueFunc pi0 scoreF isDecoyF data 
        Array.map (scoreF >> f) data 
    
    /// Calculates Q-Values from pep-values
    let getQValuesFromPEPValues (pepValues : float []) = 
        let q : float [] = Array.zeroCreate pepValues.Length
        let pep = Array.sort pepValues
        let rec loop sum i =
            if i = pepValues.Length then 
                q
            else 
                let newSum  = pep.[i] + sum
                q.[i] <- newSum / (float (i+1))
                loop newSum (i+1)
        loop 0. 0