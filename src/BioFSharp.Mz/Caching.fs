namespace BioFSharp.Mz

open System
open System.Collections.Generic

module Cache = 

    /// Wraps a SortedList to be used as a in memory cache. Especially usefull if computed theoretical fragmentation pattern
    /// or spectra are to be reused
    type Cache<[<EqualityConditionalOn; ComparisonConditionalOn >]'a,'b> = SortedList<'a,'b>

    /// Creates cache with default constructor 
    let createCache<'a,'b> = new Cache<'a,'b>()

    /// Creates cache with defined capacity of items the list can contain, gets automatically increased if the limit is reached.
    let createCacheWithCap<'a,'b> (capacity:int) = new Cache<'a,'b>(capacity)

    
    /// Creates cache with defined comparer
    let createCacheWithComp<'a,'b> (comparer:IComparer<'a>) = new Cache<'a,'b>(comparer)

    /// Creates cache with defined comparer and defined capacity of items the list can contain, gets automatically increased if the limit is reached.
    let createCacheWith<'a,'b> (comparer:IComparer<'a>) (capacity:int) = new Cache<'a,'b>(capacity, comparer)

    ///
    type Border = 
        | Upper = 1
        | Lower = 2

    ///
    let binarySearch border (cache: Cache<'a,'b>) (value:'a) = 
        let keys = cache.Keys
        let comparer = cache.Comparer
        if keys.Count = 0 then
            0
        else
            let rec loop lower upper = 
                if lower >= upper then ~~~ lower 
                else
                    let middle = lower + ((upper - lower) / 2)
                    let comparisonResult = comparer.Compare(value, keys.[middle])
                    if (comparisonResult) = 0 then
                        middle
                    elif comparisonResult < 0 then
                        loop lower (middle - 1) 
                    else
                        loop (middle + 1) upper 
            match loop 0 (keys.Count-1) with 
            | idx when idx >= 0 -> 
                idx  
            | idx when idx <  0 -> 
                let idxPos = ~~~ idx
                if idxPos >= 0 && idxPos < (keys.Count-1) then
                    if idxPos > 0 then 
                        if border = Border.Lower then
                            idxPos-1
                        else 
                            idxPos
                    else
                        idxPos
                else keys.Count-1

    /// Adds item to the Cache
    let addItem (cache: Cache<'a,'b>) (item:'a*'b) =
        match cache.ContainsKey (fst item) with
        | true     -> 
        let idx = 
            cache.IndexOfKey (fst item)
        let tempVal = cache.Values.[idx] 
        cache.RemoveAt idx 
        cache.Add item 
        | false    ->  
            cache.Add item   

    /// Adds list of items to the cache
    let bulkInsertBy (cache: Cache<'a,'b>) (items:('a*'b) list) =
        items
        |> List.map (addItem cache)

    /// Deletes list of items to the cache. Deletes all items with a key smaller than the defined cutOff.
    let bulkDeleteBy (cache: Cache<'a,'b>) (cutOff:'a) =
        let startIdx = binarySearch Border.Upper cache cutOff
        let cleanedCache = createCacheWith cache.Comparer cache.Capacity
        for i = startIdx to cache.Count-1 do
            cache.Add (cache.Keys.[i], cache.Values.[i])
        cleanedCache

    /// Returns with defined key
    let getItemBy (cache: Cache<'a,'b>) (key: 'a) = 
        cache.TryGetValue(key)
               
    /// Returns list of items with keys in the defined range. 
    let getItemsByRange (cache: Cache<'a,'b>) (range: 'a*'a) = 
        let lowerIdx = binarySearch Border.Upper cache (fst range)
        let upperIdx = binarySearch Border.Upper cache (snd range)
        [
        for i = lowerIdx to upperIdx do 
            yield (cache.Keys.[i], cache.Values.[i])
        ]

    /// Returns list of items with indices in the defined range. 
    let getItemsByIdx (cache: Cache<'a,'b>) (idxs: int*int) = 
        let lowerIdx = fst idxs
        let upperIdx = snd idxs
        [
        for i = lowerIdx to upperIdx do 
            yield (cache.Keys.[i], cache.Values.[i])
        ]

    /// Returns list of values of items with keys in the defined range. 
    let getValuesBy (cache: Cache<'a,'b>) (range: 'a*'a) = 
        let lowerIdx = binarySearch Border.Upper cache (fst range)
        let upperIdx = binarySearch Border.Upper cache (snd range)
        [        
        for i = lowerIdx to upperIdx do 
            yield (cache.Values.[i])
        ]
    
    /// Returns list of values of items with indices in the defined range. 
    let getValuesByIdx (cache: Cache<'a,'b>) (idxs: int*int) = 
        let lowerIdx = fst idxs
        let upperIdx = snd idxs
        [
        for i = lowerIdx to upperIdx do 
            yield (cache.Values.[i])
        ]
                              
    /// Returns the indices of the items in the sorted cache including items at the border of the input range of keys. 
    let containsItemsBetween (cache: Cache<'a,'b>) (range: 'a*'a) = 
        let lower,upper = range
        let keys = cache.Keys
        if keys.Count = 0 then None 
        else
        let lowerIdx = 
            let idx = binarySearch Border.Upper cache lower
            if keys.[idx] < lower then idx + 1 else idx
        let upperIdx = 
            let idx = binarySearch Border.Upper cache upper
            if keys.[idx] > upper then idx - 1 else idx
        match lowerIdx < upperIdx && lowerIdx <> upperIdx with
        | true  -> Some (lowerIdx,upperIdx)
        | false -> None    


         