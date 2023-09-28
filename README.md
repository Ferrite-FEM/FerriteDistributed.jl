# FerriteDistributed.jl

![Build Status](https://github.com/Ferrite-FEM/FerriteDistributed.jl/actions/workflows/CI.yml/badge.svg?branch=main)

This package contains the distributed assembly infrastructure for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

> [!IMPORTANT]
> This package is still experimental. Breaking changes on the design can be expected. Feedback is appreciated.

> [!NOTE]
> The current implementation is tested with [Ferrite#7e8a571](https://github.com/Ferrite-FEM/Ferrite.jl/commit/7e8a57150738094bb951d6e672fdeed205d0a1ff). This package is not compatible with the latest Ferrite release (0.3.14).

## Documentation

[![][docs-dev-img]][docs-dev-url]

## Installation
You can install FerriteDistributed from the Pkg REPL:
```
pkg> add Ferrite#7e8a571, FerriteDistributed
```

> [!NOTE]
> You need a properly configured MPI.jl installation for FerriteDistributed to work. You can consult
> the [latest MPI.jl configuration docs][mpi-config-docs] for details.

## Contributing

Contributions in all forms (bug reports, documentation, features, suggestions, ...) are very
welcome. See [CONTRIBUTING](CONTRIBUTING.md) for more details.

## Questions

If you have questions about Ferrite.jl you're welcome to reach out to us on the [Julia
Slack][julia-slack] under `#ferrite-fem` or on [Zulip][julia-zulip] under `#Ferrite.jl`.
Alternatively you can start a [new discussion][gh-discussion] in the discussion forum on the
repository. Feel free to ask us even if you are not sure the problem is with Ferrite.jl.

If you encounter what you think is a bug please report it, see
[CONTRIBUTING.md](CONTRIBUTING.md#reporting-issues) for more information.

## Community Standards

Please keep in mind that we are part of the Julia community and adhere to the
[Julia Community Standards][standards].




[docs-dev-img]: https://img.shields.io/badge/docs-dev%20release-blue
[docs-dev-url]: http://ferrite-fem.github.io/FerriteDistributed.jl/

[mpi-config-docs]: https://juliaparallel.org/MPI.jl/latest/configuration/
[standards]: https://julialang.org/community/standards/
[julia-slack]: https://julialang.org/slack/
[julia-zulip]: https://julialang.zulipchat.com/
[gh-discussion]: https://github.com/Ferrite-FEM/FerriteDistributed.jl/discussions/new
