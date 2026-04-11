# MyModule.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/MyModule.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/MyModule.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.12+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

<a id="badge-top"></a>
[![codecov](https://codecov.io/gh/sotashimozono/template.jl/graph/badge.svg?token=Q3oEEiz9A2)](https://codecov.io/gh/sotashimozono/template.jl)
[![Build Status](https://github.com/sotashimozono/MyModule.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/MyModule.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

this repository is made for template folder for developing julia project.  
some of convenient features are available, but you need to fix to your current calculations.

## TODO LIST

1. **GitHub Repository Settings**
    * [ ] **Actions Permissions**: Go to `Settings > Actions > General` and change **Workflow permissions** to **"Read and write permissions"**. This is required for `Documenter.jl` (docs) and `TagBot` to function.
    * [ ] **Allow Auto-merge**: (Recommended) Enable **"Allow auto-merge"** in `Settings > General` to streamline the PR process.
2. **Testing & Code Quality**
    * [ ] **Codecov Setup**: 
        1. Register your repository at [Codecov](https://app.codecov.io/) to obtain an upload token.
        2. Add the token to `Settings > Secrets and variables > Actions` as a repository secret named `CODECOV_TOKEN`.
        3. Replace the **[Codecov Badge](#badge-top)** link at the top of this README with the one provided in your Codecov dashboard.
3. **Documentation (Optional)**
    * [ ] **Enable Workflow**: Rename `.github/workflows/Documentation.yml.disabled` to `.github/workflows/Documentation.yml` to enable automatic document building.
    * [ ] **GitHub Pages**: After the first successful documentation build, go to `Settings > Pages` and set the source to the `gh-pages` branch.
4. **Personalization**
    * [ ] **LICENSE**: Update the year and name in the `LICENSE` file.
    * [ ] **Badges**: Ensure all badge URLs at the top of this README point to your new repository path instead of the template.

## IF YOU HAD SOME TROUBLES PLEASE MAKE `ISSUES` [HERE](https://github.com/sotashimozono/template.jl/issues)
