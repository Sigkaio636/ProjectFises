# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_A         = 200,
    N_B         = 10,
    N_C         = 5,
    Lx          = 120.0,
    Ly          = 120.0,
    T_init      = 10.0,
    dt          = 0.2,
    n_steps     = 2000,
    p_react     = 0.05,
    save_every  = 5,
    out_dir     = "reactive_gas_output",
    save_frames = false,
)