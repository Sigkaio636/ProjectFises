# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_H2O       = 600,
    N_H3O       = 600,
    N_OH        = 600,
    Lx          = 220.0,
    Ly          = 220.0,
    T_init      = 5.0,
    dt          = 0.05,
    n_steps     = 5000,
    p_react_fw  = 1.0,
    p_react_rv  = 0.5,
    save_every  = Inf, # Inf : to not generate the .gif -> just Statistics, VERY QUICK
    out_dir     = "reactive_gas_output",
    save_frames = false,
)
# EXECUTE with :>  julia --threads auto main.jl  
