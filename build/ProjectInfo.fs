module ProjectInfo

open Fake.Core

let project = "BioFSharp.Mz"

let testProject = "tests/BioFSharp.Mz.Tests/BioFSharp.Mz.Tests.fsproj"

let summary = "BioFSharp.Mz - modular computational proteomics."

let solutionFile  = "BioFSharp.Mz.sln"

let configuration = "Release"

// Git configuration (used for publishing documentation in gh-pages branch)
// The profile where the project is posted
let gitOwner = "CSBiology"

let gitName = "BioFSharp.Mz"

let gitHome = sprintf "%s/%s" "https://github.com" gitOwner

let projectRepo = sprintf "%s/%s/%s" "https://github.com" gitOwner gitName

let website = "/BioFSharp.Mz"

let pkgDir = "pkg"

/// Reads release notes that may start with an "Unreleased" section without a version. That
/// section is skipped, so the version comes from the uppermost entry that has a version.
let loadReleaseNotes (path: string) =
    let isHeading (line: string) = line.TrimStart().StartsWith "#"
    let isUnreleased (line: string) =
        isHeading line && System.String.Equals(line.Trim().TrimStart('#').Trim(), "Unreleased", System.StringComparison.OrdinalIgnoreCase)
    let lines = System.IO.File.ReadAllLines path |> Array.toList |> List.skipWhile System.String.IsNullOrWhiteSpace
    let released =
        match lines with
        | first :: rest when isUnreleased first -> rest |> List.skipWhile (isHeading >> not)
        | _ -> lines
    if List.isEmpty released then
        failwithf "%s has no entry with a version below its Unreleased section." path
    ReleaseNotes.parse released

let release = loadReleaseNotes "RELEASE_NOTES.md"

let stableVersion = SemVer.parse release.NugetVersion

let stableVersionTag = (sprintf "%i.%i.%i" stableVersion.Major stableVersion.Minor stableVersion.Patch )

let mutable prereleaseSuffix = ""

let mutable prereleaseTag = ""

let mutable isPrerelease = false
