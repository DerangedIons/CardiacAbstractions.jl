using CardiacAbstractions
using Documenter

DocMeta.setdocmeta!(
    CardiacAbstractions, :DocTestSetup, :(using CardiacAbstractions); recursive = true
)

makedocs(;
    modules = [CardiacAbstractions],
    authors = "Kyle Beggs",
    sitename = "CardiacAbstractions.jl",
    repo = Documenter.Remotes.GitHub("DerangedIons", "CardiacAbstractions.jl"),
    format = Documenter.HTML(;
        canonical = "https://DerangedIons.github.io/CardiacAbstractions.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Cell-Model Contract" => "cell_models.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/DerangedIons/CardiacAbstractions.jl",
    devbranch = "main",
    push_preview = true,
)
