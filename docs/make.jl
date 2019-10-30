using Documenter, Pulses

makedocs(;
    modules=[Pulses],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/scjethro/Pulses.jl/blob/{commit}{path}#L{line}",
    sitename="Pulses.jl",
    authors="Alastair Marshall",
    assets=String[],
)

deploydocs(;
    repo="github.com/scjethro/Pulses.jl",
)
