namespace BioFSharp.Mz

open FSharp.Stats

type PeakArray<[<EqualityConditionalOn; ComparisonConditionalOn >]'a when 'a :> IPeak> = array<'a>

module PeakArray =

    /// Builds a new PeakArray whose elements are the results of applying the given function to each 
    /// of the elements of the PeakArray.
    let map f pkarr : PeakArray<_> = 
        Array.map f pkarr

    /// Iterates the mz and intensity array and creates a Peak(mz,intensity) for each value pair. 
    /// Returns a new Peak array. 
    let zip (mz:array<float>) (intensity:array<float>) : PeakArray<_> = 
        Array.map2 (fun m i -> Peak(m,i)) mz intensity

    /// Iterates the PeakArray and unzips the value fields of each peak into two seperate arrays, the first containing 
    /// the mz values, the second the intensities. 
    let unzip (pkarr : PeakArray<_>) = 
        let n = pkarr.Length
        let mz     = Array.zeroCreate n
        let intens = Array.zeroCreate n
        for i=0 to n do
            mz.[i]     <- pkarr.[i]
            intens.[i] <- pkarr.[i]
        mz,intens
    
    /// Creates a Peak(mz,intensity) for each value pair in mzIntensity. 
    /// Returns a new PeakArray. 
    let zipMzInt (mzIntensity:array<float*float>) : PeakArray<_> = 
        Array.map (fun (m,i) -> Peak(m,i)) mzIntensity

    /// Bins peaks to their next upperIntegerMass bin
    let binToUpperIntergerMass (pkarr:PeakArray<_>) (minMassBoarder:int) (maxMassBoarder:int) = 
        let maxIndex = maxMassBoarder - minMassBoarder + 1
        let array = Array.zeroCreate (maxIndex-1)
        pkarr 
        |> Array.iter (fun p ->
            let index = int(ceil p.Mz) - minMassBoarder
            if index < maxIndex-1 && index > -1 then
                array.[index] <- max array.[index] p.Intensity)
        array

    /// Bins peaks to their nearest 1 Da bin
    let peaksToNearestUnitDaltonBin (pkarr:PeakArray<_>) (minMassBoarder:int) (maxMassBoarder:int) = 
        let maxIndex = maxMassBoarder - minMassBoarder + 1        
        let array = Array.zeroCreate (maxIndex-1)
        pkarr 
        |> Array.iter (fun p ->  
            let index = int(System.Math.Round p.Mz) - minMassBoarder
            if index < maxIndex-1 && index > -1 then
                array.[index] <- max array.[index] p.Intensity)
        array

    /// Bins peaks to their nearest 1 Da bin. Filters out peaks where the mz < minMassBoarder & > maxMassBoarder
    let peaksToNearestUnitDaltonBinVector (pkarr:PeakArray<_>) (minMassBoarder:float) (maxMassBoarder:float) =
        let minMassBoarder = int minMassBoarder
        let maxMassBoarder = int maxMassBoarder
        let maxIndex = maxMassBoarder - minMassBoarder + 1        
        let vector = Vector.create (maxIndex-1) 0.
        pkarr 
        |> Array.iter (fun p ->  
            let index = int(System.Math.Round p.Mz) - minMassBoarder
            if index < maxIndex-1 && index > -1 then
                vector.[index] <- max vector.[index] p.Intensity)
        vector
         
       