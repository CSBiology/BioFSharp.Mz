namespace BioFSharp.Mz

type PeakList<[<EqualityConditionalOn; ComparisonConditionalOn >]'a when 'a :> IPeak> = list<'a>

module PeakList =
    
    /// Builds a new PeakArray whose elements are the results of applying the given function to each 
    /// of the elements of the PeakArray.
    let map f pkl : PeakList<_> = 
        List.map f pkl

    /// Iterates the mz and intensity lists and creates a Peak(mz,intensity) for each value pair. 
    /// Returns a new PeakList. 
    let zipMzInt (mz:list<float>) (intensity:list<float>) : PeakList<_> = 
        List.map2 (fun m i -> Peak(m,i)) mz intensity

    /// Iterates the PeakList and unzips the fields of each peak into two seperate lists, the first containing 
    /// the mz values, the second the intensities. 
    let unzipMzInt (pkl : PeakList<_>) = 
        let rec loop (l:PeakList<_>) mz intens =
            match l with
            | h::t -> loop t (h.Mz::mz) (h.Intensity::intens)
            | _ -> List.rev mz ,List.rev intens           // TODO: Continuation passing
        loop pkl [] []
