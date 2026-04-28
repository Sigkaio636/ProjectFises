# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_H2O       = 300,
    N_H3O       = 300,
    N_OH        = 300,
    Lx          = 220.0,
    Ly          = 220.0,
    T_init      = 1.0,
    dt          = 0.05,
    n_steps     = 5000,
    p_react     = 1.0,
    save_every  = Inf, # Inf : to not generate the .gif -> just Statistics, VERY QUICK
    out_dir     = "reactive_gas_output",
    save_frames = false,
)
# EXECUTE with :>  julia --threads auto main.jl  
