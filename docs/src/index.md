 # LightGraphsMatching.jl: matching algorithms for LightGraphs

```@meta
CurrentModule = LightGraphsMatching
DocTestSetup = quote
    using LightGraphsMatching
    import LightGraphs
    const lg = LightGraphs
end
```

This is the documentation page for `LightGraphsMatching`. 
In all documentation examples, we assume LightGraphsMatching has been imported into scope
and that LightGraphs is available with the alias `lg`:
```julia
using LightGraphsMatching
import LightGraphs
const lg = LightGraphs
```

```@autodocs
Modules = [LightGraphsMatching]
Pages = ["LightGraphsMatching.jl", "maximum_weight_matching.jl", "lp.jl", "blossomv.jl"]
Order = [:function, :type]
```

