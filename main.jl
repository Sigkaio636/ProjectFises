# ─────────────────────────────────────────────────
# main.jl
# ─────────────────────────────────────────────────

# Incluimos el archivo que contiene la definición de tu módulo
include("src/IdealGas.jl")

# Ejecutamos la función run contenida dentro de tu módulo
IdealGas.run(
    N           = 100,
    n_steps     = 1200,
    save_every  = 5,
    radius      = 3.0,
    out_dir     = "ideal_gas_output",
    save_frames = false,
)