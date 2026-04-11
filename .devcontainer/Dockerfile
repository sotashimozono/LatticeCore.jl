FROM julia:1.12

RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN julia -e 'using Pkg; Pkg.add(["IJulia", "Revise", "JuliaFormatter"]); Pkg.precompile()'

# startup.jl（Revise）
# put startup.jl to user's home directory
# ユーザーのホームディレクトリに配置します
RUN mkdir -p /root/.julia/config && \
    echo 'try \n\
    using Revise \n\
    println("Revise loaded!") \n\
catch e \n\
    @warn "Error initializing Revise" e \n\
end' > /root/.julia/config/startup.jl

WORKDIR /workspaces
