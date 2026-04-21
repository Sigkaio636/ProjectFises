# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_H2O       = 200,
    N_H3O       = 10,
    N_OH        = 5,
    Lx          = 120.0,
    Ly          = 120.0,
    T_init      = 1.0,
    dt          = 0.2,
    n_steps     = 1500,
    p_react     = 0.05,
    save_every  = 5,
    out_dir     = "reactive_gas_output",
    save_frames = false,
)