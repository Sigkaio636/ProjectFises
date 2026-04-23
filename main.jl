# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_H2O       = 600,
    N_H3O       = 10,
    N_OH        = 10,
    Lx          = 220.0,
    Ly          = 220.0,
    T_init      = 1.0,
    dt          = 0.2,
    n_steps     = 1500,
    p_react     = 1.0,
    save_every  = 5, # Inf : to not generate the .gif -> just Statistics, VERY QUICK
    out_dir     = "reactive_gas_output",
    save_frames = false,
)