# generate examples
import Literate

EXAMPLEDIR = joinpath(@__DIR__, "src", "literate")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
mkpath(GENERATEDDIR)

# Download some assets
# include("download_resources.jl")

# Run Literate on all examples
@timeit dto "Literate." for example in readdir(EXAMPLEDIR)
    if endswith(example, ".jl")
        input = abspath(joinpath(EXAMPLEDIR, example))
        name = basename(input)
        if !liveserver
            script = @timeit dto "script()" @timeit dto name Literate.script(input, GENERATEDDIR)
            code = strip(read(script, String))
        else
            code = "<< no script output when building as draft >>"
        end

        # remove "hidden" lines which are not shown in the markdown
        line_ending_symbol = occursin(code, "\r\n") ? "\r\n" : "\n"
        code_clean = join(filter(x->!endswith(x,"#hide"),split(code, r"\n|\r\n")), line_ending_symbol)
        code_clean = replace(code_clean, r"^# This file was generated .*$"m => "")
        code_clean = strip(code_clean)

        function mdpost(str)
            str = replace(str, "@__CODE__" => code_clean)
            # MPI examples cannot run in Documenter's serial build process,
            # so convert @example blocks to static julia code blocks
            # (mimicking MPI.jl's docs approach).
            str = replace(str, r"^````@example\s+\S+\s*$"m => "````julia")
            # Remove lines that only exist to hide output in @example blocks
            str = replace(str, r"^nothing\s+#\s*hide\s*\n"m => "")
            return str
        end

        @timeit dto "markdown()" @timeit dto name begin
            Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
        end
    elseif any(endswith.(example, [".png", ".jpg", ".gif"]))
        cp(joinpath(EXAMPLEDIR, example), joinpath(GENERATEDDIR, example); force=true)
    else
        @warn "ignoring $example"
    end
end

# remove any .vtu files in the generated dir (should not be deployed)
@timeit dto "remove vtk files" cd(GENERATEDDIR) do
    foreach(file -> endswith(file, ".pvtu") && rm(file), readdir())
    foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
    foreach(file -> endswith(file, ".pvd") && rm(file), readdir())
end