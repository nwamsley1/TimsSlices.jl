# A .d bundle opened for reading: SQLite tables in memory, analysis.tdf_bin memory-mapped.
using Mmap

struct TdfFile
    dir::String
    meta::Dict{String, String}
    frames::FrameTable
    dia::DiaScheme
    mz_cal::LinearMzCal
    mz_cal_resid_ppm::Vector{Float64}
    im_cal::LinearImCal
    ce_ramp::CeRamp
    bin::Vector{UInt8}          # mmap of analysis.tdf_bin
    compression::Int
    max_scans::Int
    max_peaks::Int
end

function open_tdf(dir::AbstractString)
    isdir(dir) || error("not a directory: $dir")
    db = SQLite.DB(joinpath(dir, "analysis.tdf"))
    meta = read_global_metadata(db)
    frames = read_frames(db)
    dia = read_dia_scheme(db, frames)
    mzcal_id = isempty(frames.mz_calibration) ? 1 : Int(frames.mz_calibration[1])
    cal, resid = regressed_mz_cal(db, mzcal_id)
    imcal = boundary_im_cal(meta, frames)
    ce = fit_ce_ramp(dia)
    close(db)
    compression = parse(Int, get(meta, "TimsCompressionType", "2"))
    compression == 2 || error("only TimsCompressionType 2 is supported (file has $compression)")
    bin = open(joinpath(dir, "analysis.tdf_bin"), "r") do io
        Mmap.mmap(io, Vector{UInt8}, filesize(io))
    end
    TdfFile(String(dir), meta, frames, dia, cal, resid, imcal, ce, bin, compression,
            isempty(frames.num_scans) ? 0 : Int(maximum(frames.num_scans)),
            isempty(frames.num_peaks) ? 0 : Int(maximum(frames.num_peaks)))
end

n_frames(f::TdfFile) = length(f.frames)
is_ms1(f::TdfFile, i::Integer) = f.frames.msms_type[i] == 0
is_dia(f::TdfFile, i::Integer) = f.frames.msms_type[i] == 9
"Rows of the frames that are converted: MS1 and diaPASEF MS2 frames, in order."
valid_frames(f::TdfFile) = [i for i in 1:n_frames(f) if is_ms1(f, i) || is_dia(f, i)]
"Windows of frame row i (MS1: one window over all scans)."
function windows(f::TdfFile, i::Integer)
    if is_ms1(f, i)
        return (DiaWindow(0, f.frames.num_scans[i], 0f0, 0f0, 0f0),)
    else
        g = f.dia.frame_group[i]
        g == 0 && error("frame $(f.frames.id[i]) has MsMsType 9 but no DiaFrameMsMsInfo row")
        return f.dia.groups[g]
    end
end
"Number of TOF bins on the digitiser axis (upper bound for any tof value)."
n_bins(f::TdfFile) = parse(Int, get(f.meta, "DigitizerNumSamples", "0"))

"Raw block header at frame row i: (payload view, block_size, scan_count)."
function raw_block(f::TdfFile, i::Integer)
    off = f.frames.tims_id[i]
    off + 8 <= length(f.bin) || error("frame $(f.frames.id[i]): block offset $off beyond file")
    p = pointer(f.bin) + off
    block_size = GC.@preserve f unsafe_load(Ptr{UInt32}(p))
    scan_count = GC.@preserve f unsafe_load(Ptr{UInt32}(p + 4))
    off + block_size <= length(f.bin) || error("frame $(f.frames.id[i]): block runs beyond file")
    view(f.bin, off + 9:off + block_size), Int(block_size), Int(scan_count)
end

"Decode frame row i into `buf`."
function read_frame!(buf::FrameBuffer, f::TdfFile, i::Integer)
    payload, _, scan_count = raw_block(f, i)
    ns = Int(f.frames.num_scans[i]); np = Int(f.frames.num_peaks[i])
    scan_count == ns || error("frame $(f.frames.id[i]): block scan_count $scan_count != NumScans $ns")
    decode_codec2!(buf, payload, ns, np)
end
