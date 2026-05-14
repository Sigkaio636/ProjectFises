# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Carga del módulo principal
include("src/IdealGas.jl")

# Ejecución de la simulación con parámetros personalizados
IdealGas.run(
    N_H2O           = 300,
    N_H3O           = 600,
    N_OH            = 1,
    Lx              = 200.0,
    Ly              = 200.0,
    T_init          = 1.0,
    dt              = 0.05,
    n_steps         = 8000,
    p_react_fw      = 1.0,
    p_react_rv      = 0.5,
    add_spe_every   = 100,
    amount_spe      = (0,0,2),
    save_every      = Inf, # Inf : to not generate the .gif -> just Statistics, VERY QUICK
    out_dir         = "valoration_output",
    save_frames     = false,
)
# EXECUTE with :>  julia --threads auto main.jl  
