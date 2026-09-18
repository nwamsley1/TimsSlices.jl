using Test, Random, Statistics
using TimsSlices
const TS = TimsSlices

# Real-data tests need the PRIDE files; set TIMSSLICES_TEST_DATA to the directory holding the .d bundles
# (default ~/BrukerTims/pride). TIMSSLICES_TEST_BIG=1 also runs the equivalence test on the E. coli and
# 250 pg files (minutes).
const DATA = get(ENV, "TIMSSLICES_TEST_DATA", expanduser("~/BrukerTims/pride"))
const HELA = joinpath(DATA, "20210510_TIMS03_EVO03_PaSk_SA_HeLa_50ng_5_6min_DIA_high_speed_S1-B2_1_25186.d")
const ECOLI = joinpath(DATA, "LFQ_Ultra2_diaPASEF_5min_50ng_Ecoli_01.d")
const PG250 = joinpath(DATA, "LFQ_Ultra_diaPASEF_15min_250pg_Human_01.d")
const BIG = get(ENV, "TIMSSLICES_TEST_BIG", "0") == "1"

@testset "TimsSlices" begin
    include("test_codec.jl")
    include("test_smooth.jl")
    if isdir(HELA)
        include("test_realdata.jl")
    else
        @warn "HeLa test file not found at $HELA; skipping real-data tests"
    end
end
