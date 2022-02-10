using pylib::FileReadString
using std::fs::File
using std::fs::OpenOptions
import astropy.constants
import astropy.units
import pandas
import numpy

using scipy: optimize
import PyMieScatt
using root_functions: advdiff, vfall, vfall_find_root, qvs_below_model, find_cond_t, solve_force_balance
using calc_mie: fort_mie_calc, calc_new_mieff
using ::: gas_properties
using ::: pvaps
using justplotit: plot_format, find_nearest_1d
using direct_mmr_solver: direct_solver
function compute{T0, T1, T2, T3, T4, T5, T6, T7, T8, RT}(atmo::T0, directory::T1, as_dict::T2, og_solver::T3, direct_tol::T4, refine_TP::T5, og_vfall::T6, analytical_rg::T7, do_virtual::T8)::RT
mmw = atmo.mmw
mh = atmo.mh
condensibles = atmo.condensibles
ngas = length(condensibles)
gas_mw = zeros(np, ngas)
gas_mmr = zeros(np, ngas)
rho_p = zeros(np, ngas)
H = ((atmo.r_atmos*atmo.Teff)/atmo.g)
for (i, igas) in zip((0:ngas - 1), condensibles)
run_gas = getattr(gas_properties, igas)
gas_mw[i], gas_mmr[i], rho_p[i] = run_gas(mmw, mh, atmo.gas_mmr)
qext_gas, qscat_gas, cos_qscat_gas, nwave, radius, wave_in = get_mie(igas, directory)
if i == 0
nradii = length(radius)
rmin = min(np, radius)
radius, rup, dr = get_r_grid(rmin)
qext = zeros(np, (nwave, nradii, ngas))
qscat = zeros(np, (nwave, nradii, ngas))
cos_qscat = zeros(np, (nwave, nradii, ngas))
end
qext[None], qscat[None], cos_qscat[None] = (qext_gas, qscat_gas, cos_qscat_gas)
end
if og_solver
if atmo.param == "exp"
atmo.b = ((6*atmo.b)*H)
fsed_in = (atmo.fsed - atmo.eps)
else

if atmo.param == "const"
fsed_in = atmo.fsed;
end
end
qc, qt, rg, reff, ndz, qc_path, mixl, z_cld = eddysed(atmo.t_level, atmo.p_level, atmo.t_layer, atmo.p_layer, condensibles, gas_mw, gas_mmr, rho_p, mmw, atmo.g, atmo.kz, atmo.mixl, fsed_in, atmo.b, atmo.eps, atmo.scale_h, atmo.z_top, atmo.z_alpha, min(atmo.z), atmo.param, mh, atmo.sig, rmin, nradii, atmo.d_molecule, atmo.eps_k, atmo.c_p_factor, og_vfall)
pres_out = atmo.p_layer
temp_out = atmo.t_layer
z_out = atmo.z
else

fsed_in = atmo.fsed
z_cld = nothing
qc, qt, rg, reff, ndz, qc_path, pres_out, temp_out, z_out, mixl = direct_solver(atmo.t_layer, atmo.p_layer, condensibles, gas_mw, gas_mmr, rho_p, mmw, atmo.g, atmo.kz, atmo.fsed, mh, atmo.sig, rmin, nradii, atmo.d_molecule, atmo.eps_k, atmo.c_p_factor, direct_tol, refine_TP, og_vfall, analytical_rg)
end
opd, w0, g0, opd_gas = calc_optics(nwave, qc, qt, rg, reff, ndz, radius, dr, qext, qscat, cos_qscat, atmo.sig, rmin, nradii)
if as_dict
if atmo.param == "exp"
fsed_out = ((fsed_in*exp(np, ((atmo.z - atmo.z_alpha)/atmo.b))) + atmo.eps)
else

fsed_out = fsed_in
end
return create_dict(qc, qt, rg, reff, ndz, opd, w0, g0, opd_gas, wave_in, pres_out, temp_out, condensibles, mh, mmw, fsed_out, atmo.sig, nradii, rmin, z_out, atmo.dz_layer, mixl, atmo.kz, atmo.scale_h, z_cld)
else

return (opd, w0, g0)
end
end

function create_dict{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, RT}(qc::T0, qt::T1, rg::T2, reff::T3, ndz::T4, opd::T5, w0::T6, g0::T7, opd_gas::T8, wave::T9, pressure::T10, temperature::T11, gas_names::T12, mh::T13, mmw::T14, fsed::T15, sig::T16, nrad::T17, rmin::T18, z::T19, dz_layer::T20, mixl::T21, kz::T22, scale_h::T23, z_cld::T24)::RT
return Dict("pressure" => (pressure/1000000.0), "pressure_unit" => "bar", "temperature" => temperature, "temperature_unit" => "kelvin", "wave" => wave[None], "wave_unit" => "micron", "condensate_mmr" => qc, "cond_plus_gas_mmr" => qt, "mean_particle_r" => (rg*10000.0), "droplet_eff_r" => (reff*10000.0), "r_units" => "micron", "column_density" => ndz, "column_density_unit" => "#/cm^2", "opd_per_layer" => opd, "single_scattering" => w0, "asymmetry" => g0, "opd_by_gas" => opd_gas, "condensibles" => gas_names, "scalar_inputs" => Dict("mh" => mh, "mmw" => mmw, "sig" => sig, "nrad" => nrad, "rmin" => rmin), "fsed" => fsed, "altitude" => z, "layer_thickness" => dz_layer, "z_unit" => "cm", "mixing_length" => mixl, "mixing_length_unit" => "cm", "kz" => kz, "kz_unit" => "cm^2/s", "scale_height" => scale_h, "cloud_deck" => z_cld)
end

function calc_optics{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, RT}(nwave::T0, qc::T1, qt::T2, rg::T3, reff::T4, ndz::T5, radius::T6, dr::T7, qext::T8, qscat::T9, cos_qscat::T10, sig::T11, rmin::T12, nrad::T13, verbose::T14)::RT
PI = np.pi
nz = qc.shape[0]
ngas = qc.shape[1]
nrad = length(radius);
opd_layer = zeros(np, (nz, ngas))
scat_gas = zeros(np, (nz, nwave, ngas))
ext_gas = zeros(np, (nz, nwave, ngas))
cqs_gas = zeros(np, (nz, nwave, ngas))
opd = zeros(np, (nz, nwave))
opd_gas = zeros(np, (nz, ngas))
w0 = zeros(np, (nz, nwave))
g0 = zeros(np, (nz, nwave))
warning = ""
for iz in (0:nz - 1)
for igas in (0:ngas - 1)
if ndz[(iz, igas)] > 0
if log10(np, rg[(iz, igas)]) < (log10(np, rmin) + (0.75*sig))
warning0 = join("", ["Take caution in analyzing results. There have been a calculated particle radii off the Mie grid, which has a min radius of ", string(rmin), "cm and distribution of ", string(sig), ". The following errors:"])
warning += format("{0}cm for the {1}th gas at the {2}th grid point; ", string(rg[(iz, igas)]), string(igas), string(iz))
end
r2 = (pow(rg[(iz, igas)], 2)*exp(np, (2*pow(np.log(sig), 2))))
opd_layer[(iz, igas)] = (((2.0*PI)*r2)*ndz[(iz, igas)])
rsig = sig
norm = 0.0
for irad in (0:nrad - 1)
rr = radius[irad]
arg1 = (dr[irad]/((sqrt(np, (2.0*PI))*rr)*log(np, rsig)))
arg2 = (-(pow(log(np, (rr/rg[(iz, igas)])), 2))/(2*pow(log(np, rsig), 2)))
norm = (norm + (arg1*exp(np, arg2)));
end
norm = (ndz[(iz, igas)]/norm);
for irad in (0:nrad - 1)
rr = radius[irad]
arg1 = (dr[irad]/(sqrt(np, (2.0*PI))*log(np, rsig)))
arg2 = (-(pow(log(np, (rr/rg[(iz, igas)])), 2))/(2*pow(log(np, rsig), 2)))
pir2ndz = ((((norm*PI)*rr)*arg1)*exp(np, arg2))
for iwave in (0:nwave - 1)
scat_gas[(iz, iwave, igas)] = (scat_gas[(iz, iwave, igas)] + (qscat[(iwave, irad, igas)]*pir2ndz))
ext_gas[(iz, iwave, igas)] = (ext_gas[(iz, iwave, igas)] + (qext[(iwave, irad, igas)]*pir2ndz))
cqs_gas[(iz, iwave, igas)] = (cqs_gas[(iz, iwave, igas)] + (cos_qscat[(iwave, irad, igas)]*pir2ndz))
end
end
end
end
end
for iz in (0:nz - 1)
for iwave in (0:nwave - 1)
opd_scat = 0.0
opd_ext = 0.0
cos_qs = 0.0
for igas in (0:ngas - 1)
opd_scat = (opd_scat + scat_gas[(iz, iwave, igas)]);
opd_ext = (opd_ext + ext_gas[(iz, iwave, igas)]);
cos_qs = (cos_qs + cqs_gas[(iz, iwave, igas)]);
if opd_scat > 0.0
opd[(iz, iwave)] = opd_ext
w0[(iz, iwave)] = (opd_scat/opd_ext)
g0[(iz, iwave)] = (cos_qs/opd_scat)
end
end
end
end
opd_tot = 0.0
for igas in (0:ngas - 1)
opd_gas[(0, igas)] = opd_layer[(0, igas)]
for iz in (1:nz - 1)
opd_gas[(iz, igas)] = (opd_gas[((iz - 1), igas)] + opd_layer[(iz, igas)])
end
end
if (warning != "" & verbose)
println(join([((warning0 + warning) + " Turn off warnings by setting verbose=False.")], " "));
end
return (opd, w0, g0, opd_gas)
end

function eddysed{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, RT}(t_top::T0, p_top::T1, t_mid::T2, p_mid::T3, condensibles::T4, gas_mw::T5, gas_mmr::T6, rho_p::T7, mw_atmos::T8, gravity::T9, kz::T10, mixl::T11, fsed::T12, b::T13, eps::T14, scale_h::T15, z_top::T16, z_alpha::T17, z_min::T18, param::T19, mh::T20, sig::T21, rmin::T22, nrad::T23, d_molecule::T24, eps_k::T25, c_p_factor::T26, og_vfall::T27, do_virtual::T28, supsat::T29, verbose::T30)::RT
did_gas_condense = condensibles.iter().map(|i| false).collect::<Vec<_>>()
t_bot = t_top[-1]
p_bot = p_top[-1]
z_bot = z_top[-1]
ngas = length(condensibles)
nz = length(t_mid)
qc = zeros(np, (nz, ngas))
qt = zeros(np, (nz, ngas))
rg = zeros(np, (nz, ngas))
reff = zeros(np, (nz, ngas))
ndz = zeros(np, (nz, ngas))
fsed_layer = zeros(np, (nz, ngas))
qc_path = zeros(np, ngas)
z_cld_out = zeros(np, ngas)
for (i, igas) in zip((0:ngas - 1), condensibles)
q_below = gas_mmr[i]
if do_virtual
z_cld = nothing
qvs_factor = (((supsat + 1)*gas_mw[i])/mw_atmos)
get_pvap = getattr(pvaps, igas)
if igas == "Mg2SiO4"
pvap = get_pvap(t_bot, p_bot, mh)
else

pvap = get_pvap(t_bot, mh)
end
qvs = ((qvs_factor*pvap)/p_bot)
if qvs <= q_below
p_lo = p_bot
p_hi = (p_bot*1000.0)
dtdlnp = ((t_top[-2] - t_bot)/log(np, (p_bot/p_top[-2])))
qv_dtdlnp = dtdlnp
qv_p = p_bot
qv_t = t_bot
qv_gas_name = igas
qv_factor = qvs_factor
" try:

                    p_base = optimize.root_scalar(qvs_below_model, 
                                bracket=[p_lo, p_hi], method='brentq', 
                                args=(qv_dtdlnp,qv_p, qv_t,qv_factor ,qv_gas_name,mh,q_below)
                                )#, xtol = 1e-20)

                    if verbose: print('Virtual Cloud Found: '+ qv_gas_name)
                    root_was_found = True
                except ValueError: 
                    root_was_found = False
                ";
if root_was_found
did_gas_condense[i] = true
p_base = p_base.root
t_base = (t_bot + (log(np, (p_bot/p_base))*dtdlnp))
z_base = (z_bot + (scale_h[-1]*log(np, (p_bot/p_base))))
p_layer_virtual = (0.5*(p_bot + p_base))
t_layer_virtual = (t_bot + (log10(np, (p_bot/p_layer_virtual))*dtdlnp))
qc_v, qt_v, rg_v, reff_v, ndz_v, q_below, z_cld, fsed_layer_v = layer(igas, rho_p[i], t_layer_virtual, convert(, p_layer_virtual), t_bot, t_base, p_bot, p_base, kz[-1], mixl[-1], gravity, mw_atmos, gas_mw[i], q_below, supsat, fsed, b, eps, z_bot, z_base, z_alpha, z_min, param, sig, mh, rmin, nrad, d_molecule, eps_k, c_p_factor, og_vfall, z_cld)
end
end
end
z_cld = nothing
for iz in ((nz - 1):-1:-1-1)
qc[(iz, i)], qt[(iz, i)], rg[(iz, i)], reff[(iz, i)], ndz[(iz, i)], q_below, z_cld, fsed_layer[(iz, i)] = layer(igas, rho_p[i], t_mid[iz], p_mid[iz], t_top[iz], t_top[(iz + 1)], p_top[iz], p_top[(iz + 1)], kz[iz], mixl[iz], gravity, mw_atmos, gas_mw[i], q_below, supsat, fsed, b, eps, z_top[iz], z_top[(iz + 1)], z_alpha, z_min, param, sig, mh, rmin, nrad, d_molecule, eps_k, c_p_factor, og_vfall, z_cld)
qc_path[i] = (qc_path[i] + ((qc[(iz, i)]*(p_top[(iz + 1)] - p_top[iz]))/gravity))
end
z_cld_out[i] = z_cld
end
return (qc, qt, rg, reff, ndz, qc_path, mixl, z_cld_out)
end

function layer{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, RT}(gas_name::T0, rho_p::T1, t_layer::T2, p_layer::T3, t_top::T4, t_bot::T5, p_top::T6, p_bot::T7, kz::T8, mixl::T9, gravity::T10, mw_atmos::T11, gas_mw::T12, q_below::T13, supsat::T14, fsed::T15, b::T16, eps::T17, z_top::T18, z_bot::T19, z_alpha::T20, z_min::T21, param::T22, sig::T23, mh::T24, rmin::T25, nrad::T26, d_molecule::T27, eps_k::T28, c_p_factor::T29, og_vfall::T30, z_cld::T31)::RT
nsub_max = 128
R_GAS = 83143000.0
AVOGADRO = 6.02e+23
K_BOLTZ = (R_GAS/AVOGADRO)
PI = np.pi
nsub = 1
r_atmos = (R_GAS/mw_atmos)
r_cloud = (R_GAS/gas_mw)
c_p = (c_p_factor*r_atmos)
dp_layer = (p_bot - p_top)
dlnp = log(np, (p_bot/p_top))
dtdlnp = ((t_top - t_bot)/dlnp)
lapse_ratio = (((t_bot - t_top)/dlnp)/(t_layer/c_p_factor))
rho_atmos = (p_layer/(r_atmos*t_layer))
scale_h = ((r_atmos*t_layer)/gravity)
w_convect = (kz/mixl)
n_atmos = (p_layer/(K_BOLTZ*t_layer))
mfp = (1.0/(((sqrt(np, 2.0)*n_atmos)*PI)*pow(d_molecule, 2)))
visc = ((((5.0/16.0)*sqrt(np, (((PI*K_BOLTZ)*t_layer)*(mw_atmos/AVOGADRO))))/(PI*pow(d_molecule, 2)))/(1.22*pow((t_layer/eps_k), -0.16)))
converge = false
while !(converge)
qc_layer = 0.0
qt_layer = 0.0
ndz_layer = 0.0
opd_layer = 0.0
qt_bot_sub = q_below
p_bot_sub = p_bot
z_bot_sub = z_bot
dp_sub = (dp_layer/nsub)
for isub in (0:nsub - 1)
qt_below = qt_bot_sub
p_top_sub = (p_bot_sub - dp_sub)
dz_sub = (scale_h*log(np, (p_bot_sub/p_top_sub)))
p_sub = (0.5*(p_bot_sub + p_top_sub))
z_top_sub = (z_bot_sub + dz_sub)
z_sub = (z_bot_sub + (scale_h*log(np, (p_bot_sub/p_sub))))
t_sub = (t_bot + (log(np, (p_bot/p_sub))*dtdlnp))
qt_top, qc_sub, qt_sub, rg_sub, reff_sub, ndz_sub, z_cld, fsed_layer = calc_qc(gas_name, supsat, t_sub, convert(, p_sub), convert(, r_atmos), convert(, r_cloud), qt_below, mixl, convert(, dz_sub), gravity, mw_atmos, convert(, mfp), convert(, visc), rho_p, w_convect, fsed, b, eps, param, convert(, z_bot_sub), convert(, z_sub), z_alpha, z_min, sig, mh, rmin, nrad, og_vfall, z_cld)
qc_layer = (qc_layer + ((qc_sub*dp_sub)/gravity));
qt_layer = (qt_layer + ((qt_sub*dp_sub)/gravity));
ndz_layer = (ndz_layer + ndz_sub);
if reff_sub > 0.0
opd_layer = (opd_layer + ((((1.5*qc_sub)*dp_sub)/gravity)/(rho_p*reff_sub)));
end
qt_bot_sub = qt_top;
p_bot_sub = p_top_sub;
z_bot_sub = z_top_sub;
end
if nsub_max == 1
converge = true;
else

if nsub == 1
opd_test = opd_layer
else

if opd_layer == 0.0||nsub >= nsub_max
converge = true;
else

if abs((1.0 - (opd_test/opd_layer))) <= 0.01
converge = true;
else

opd_test = opd_layer;
end
end
end
end
nsub = (nsub*2);
end
q_below = qt_top;
if opd_layer > 0.0
reff_layer = ((1.5*qc_layer)/(rho_p*opd_layer))
lnsig2 = (0.5*pow(log(np, sig), 2))
rg_layer = (reff_layer*exp(np, (-5*lnsig2)))
else

reff_layer = 0.0
rg_layer = 0.0
end
qc_layer = ((qc_layer*gravity)/dp_layer);
qt_layer = ((qt_layer*gravity)/dp_layer);
return (qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, q_below, z_cld, fsed_layer)
end

function calc_qc{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, RT}(gas_name::T0, supsat::T1, t_layer::T2, p_layer::T3, r_atmos::T4, r_cloud::T5, q_below::T6, mixl::T7, dz_layer::T8, gravity::T9, mw_atmos::T10, mfp::T11, visc::T12, rho_p::T13, w_convect::T14, fsed::T15, b::T16, eps::T17, param::T18, z_bot::T19, z_layer::T20, z_alpha::T21, z_min::T22, sig::T23, mh::T24, rmin::T25, nrad::T26, og_vfall::T27, z_cld::T28)::RT
get_pvap = getattr(pvaps, gas_name)
if gas_name == "Mg2SiO4"
pvap = get_pvap(t_layer, p_layer, mh)
else

pvap = get_pvap(t_layer, mh)
end
fs = (supsat + 1)
rho_atmos = (p_layer/(r_atmos*t_layer))
qvs = (((fs*pvap)/(r_cloud*t_layer))/rho_atmos)
if q_below < qvs
qt_layer = q_below
qt_top = q_below
qc_layer = 0.0
rg_layer = 0.0
reff_layer = 0.0
ndz_layer = 0.0
z_cld = z_cld;
fsed_mid = 0
else

if isinstance(z_cld, type_(nothing))
z_cld = z_bot
else

z_cld = z_cld
end
qhi = q_below
qlo = (qhi/1000.0)
ad_qbelow = q_below
ad_qvs = qvs
ad_mixl = mixl
ad_dz = dz_layer
ad_rainf = fsed
if param == "const"
qt_top = (qvs + ((q_below - qvs)*exp(np, ((-(fsed)*dz_layer)/mixl))));
else

if param == "exp"
fs = (fsed/exp(np, (z_alpha/b)));
qt_top = (qvs + ((q_below - qvs)*exp(np, (((((-(b)*fs)/mixl)*np.exp((z_bot/b)))*(np.exp((dz_layer/b)) - 1)) + ((eps*dz_layer)/mixl)))));
end
end
qt_layer = (0.5*(q_below + qt_top))
qc_layer = max(np, [0.0, (qt_layer - qvs)])
rlo = 1e-10
rhi = 10.0
find_root = true
"while find_root:
            try:
                if og_vfall:
                    rw_temp = optimize.root_scalar(vfall_find_root, bracket=[rlo, rhi], method='brentq', 
                            args=(gravity,mw_atmos,mfp,visc,t_layer,p_layer, rho_p,w_convect))
                else:
                    rw_temp = solve_force_balance(\"rw\", w_convect, gravity, mw_atmos, mfp,
                                                    visc, t_layer, p_layer, rho_p, rlo, rhi)
                find_root = False
            except ValueError:
                rlo = rlo/10
                rhi = rhi*10
        ";
if og_vfall
rw_layer = rw_temp.root
else

rw_layer = rw_temp
end
lnsig2 = (0.5*pow(log(np, sig), 2))
sig_alpha = max(np, [1.1, sig])
function pow_law{T0, T1, RT}(r::T0, alpha::T1)::RT
return (log(np, w_convect) + (alpha*log(np, (r/rw_layer))))
end

r_, rup, dr = get_r_grid()
vfall_temp = []
for j in (0:length(r_) - 1)
if og_vfall
push!(vfall_temp, vfall(r_[j], gravity, mw_atmos, mfp, visc, t_layer, p_layer, rho_p));
else

vlo = 1.0
vhi = 1000000.0
find_root = true;
while find_root
"try:
                        vfall_temp.append(solve_force_balance(\"vfall\", r_[j], gravity, mw_atmos, 
                            mfp, visc, t_layer, p_layer, rho_p, vlo, vhi))
                        find_root = False
                    except ValueError:
                        vlo = vlo/10
                        vhi = vhi*10";
end
end
end
pars, cov = curve_fit(optimize, pow_law, r_, np.log(vfall_temp), [0], (-(np.inf), np.inf))
alpha = pars[0]
if param == "exp"
fsed_mid = ((fs*exp(np, (z_layer/b))) + eps)
else

fsed_mid = fsed
end
rg_layer = ((pow(fsed_mid, (1.0/alpha))*rw_layer)*exp(np, (-((alpha + 6))*lnsig2)))
reff_layer = (rg_layer*exp(np, (5*lnsig2)))
ndz_layer = (((((3*rho_atmos)*qc_layer)*dz_layer)/(((4*np.pi)*rho_p)*pow(rg_layer, 3)))*exp(np, (-9*lnsig2)))
end
return (qt_top, qc_layer, qt_layer, rg_layer, reff_layer, ndz_layer, z_cld, fsed_mid)
end

struct Atmosphere
mh::
mmw::
condensibles::
fsed::
b::
sig::
param::
eps::
verbose::
supsat::
gas_mmr::
eps_k::Float64
d_molecule::Float64
c_p_factor::Float64
R_GAS::Float64
AVOGADRO::Float64
K_BOLTZ::
p_level::Float64
t_level::
alpha_pressure::
Teff::
r_atmos::
p_layer::Float64
dtdlnp::
t_layer::
lapse_ratio::
rho_atmos::
scale_h::
c_p::
dz_pmid::
dz_layer::
z_top::
z::
z_alpha::
mixl::
kz::Float64
chf::
g::
gravity_unit::String
end

function __init__{T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10}(self::Atmosphere, condensibles::T0, fsed::T1, b::T2, eps::T3, mh::T4, mmw::T5, sig::T6, param::T7, verbose::T8, supsat::T9, gas_mmr::T10)
self.mh = mh
self.mmw = mmw
self.condensibles = condensibles
self.fsed = fsed
self.b = b
self.sig = sig
self.param = param
self.eps = eps
self.verbose = verbose
constants(self);
self.supsat = supsat
self.gas_mmr = gas_mmr
end

function constants(self::Atmosphere)
self.eps_k = 59.7
self.d_molecule = 2.827e-08
self.c_p_factor = (7.0/2.0)
self.R_GAS = 83143000.0
self.AVOGADRO = 6.02e+23
self.K_BOLTZ = (self.R_GAS/self.AVOGADRO)
end

function ptk{T0, T1, T2, T3, T4, T5, T6, T7}(self::Atmosphere, df::T0, filename::T1, kz_min::T2, constant_kz::T3, latent_heat::T4, convective_overshoot::T5, Teff::T6, alpha_pressure::T7)
if !(isinstance(df, type_(nothing)))
if isinstance(df, dict)
df = DataFrame(pd, df)
end
df = sort_values(df, "pressure");
else

if !(isinstance(filename, type_(nothing)))
df = read_csv(pd, filename, pd_kwargs);
df = sort_values(df, "pressure");
end
end
self.p_level = (array(np, df["pressure"])*1000000.0)
self.t_level = array(np, df["temperature"])
if alpha_pressure == nothing
self.alpha_pressure = min(df["pressure"])
else

self.alpha_pressure = alpha_pressure
end
get_atmo_parameters(self);
get_kz_mixl(self, df, constant_kz, latent_heat, convective_overshoot, kz_min);
if isinstance(Teff, type_(nothing))
onebar = argmin(np.abs(((self.p_level/1000000.0) - 1.0)))
self.Teff = self.t_level[onebar]
else

self.Teff = Teff
end
end

function get_atmo_parameters(self::Atmosphere)
self.r_atmos = (self.R_GAS/self.mmw)
dlnp = log(np, (self.p_level[1..]/self.p_level[0..-1]))
self.p_layer = (0.5*(self.p_level[1..] + self.p_level[0..-1]))
self.dtdlnp = ((self.t_level[0..-1] - self.t_level[1..])/dlnp)
self.t_layer = (self.t_level[1..] + (log(np, (self.p_level[1..]/self.p_layer))*self.dtdlnp))
self.lapse_ratio = (((self.t_level[1..] - self.t_level[0..-1])/dlnp)/(self.t_layer/self.c_p_factor))
self.rho_atmos = (self.p_layer/(self.r_atmos*self.t_layer))
self.scale_h = ((self.r_atmos*self.t_layer)/self.g)
self.c_p = (self.c_p_factor*self.r_atmos)
self.dz_pmid = (self.scale_h*log(np, (self.p_level[1..]/self.p_layer)))
self.dz_layer = (self.scale_h*dlnp)
self.z_top = concatenate(np, ([0], np.cumsum(self.dz_layer[..])))[..]
self.z = (self.z_top[1..] + self.dz_pmid)
p_alpha = find_nearest_1d((self.p_layer/1000000.0), self.alpha_pressure)
z_temp = cumsum(np, self.dz_layer[..])[..]
self.z_alpha = z_temp[p_alpha]
end

function get_kz_mixl{T0, T1, T2, T3, T4}(self::Atmosphere, df::T0, constant_kz::T1, latent_heat::T2, convective_overshoot::T3, kz_min::T4)
if latent_heat
self.mixl = (array(np, self.lapse_ratio.iter().map(|ilr| np.max([0.1, ilr])).collect::<Vec<_>>())*self.scale_h)
else

self.mixl = (1*self.scale_h)
end
if "kz" in keys(df)
if df.loc[df["kz"] < kz_min].shape[0] > 0
df.loc[df["kz"] < kz_min] = kz_min
if self.verbose
println(join(["Overwriting some Kz values to minimum value set by kz_min 
                     You can always turn off these warnings by setting verbose=False"], " "));
end
end
kz_level = array(np, df["kz"])
self.kz = (0.5*(kz_level[1..] + kz_level[0..-1]))
self.chf = nothing
else

if !(isinstance(constant_kz, type_(nothing)))
self.kz = (zeros(np, (df.shape[0] - 1)) + constant_kz)
self.chf = nothing
else

if "chf" in keys(df)
self.chf = array(np, df["chf"])
if !(isinstance(convective_overshoot, type_(nothing)))
used = false
nz = length(self.p_layer)
for iz in ((nz - 1):-1:-1-1)
ratio_min = ((convective_overshoot*self.p_level[iz])/self.p_level[(iz + 1)])
if self.chf[iz] < (ratio_min*self.chf[(iz + 1)])
self.chf[iz] = (self.chf[(iz + 1)]*ratio_min)
used = true;
end
end
if self.verbose
println(join(["Convective overshoot was turned on. The convective heat flux 
                    has been adjusted such that it is not allowed to decrease more than {0} 
                    the pressure. This number is set with the convective_overshoot parameter. 
                    It can be disabled with convective_overshoot=None. To turn
                    off these messages set verbose=False in Atmosphere".format(convective_overshoot)], " "));
end
end
gc_kzz = ((((1.0/3.0)*self.scale_h)*pow((self.mixl/self.scale_h), (4.0/3.0)))*pow(((self.r_atmos*self.chf[1..])/(self.rho_atmos*self.c_p)), (1.0/3.0)))
self.kz = gc_kzz.iter().map(|i| max(np, [i, kz_min])).collect::<Vec<_>>()
else

raise!(Exception("Users can define kz by: 
             1) Adding 'kz' as a column or key to your dataframe dict, or file 
             2) Defining constant-w-altitude kz through the constant_kz input 
              3) Adding 'chf', the conective heat flux as a column to your             dataframe, dict or file.")) # unsupported
end
end
end
end

function gravity{T0, T1, T2, T3, T4, T5}(self::Atmosphere, gravity::T0, gravity_unit::T1, radius::T2, radius_unit::T3, mass::T4, mass_unit::T5)
if mass != nothing&&radius != nothing
m = to((mass*mass_unit), u.g)
r = to((radius*radius_unit), u.cm)
g = ((c.G.cgs*m)/pow(r, 2)).value
self.g = g
self.gravity_unit = "cm/(s**2)"
else

if gravity != nothing
g = to((gravity*gravity_unit), "cm/(s**2)");
g = g.value;
self.g = g
self.gravity_unit = "cm/(s**2)"
else

raise!(Exception("Need to specify gravity or radius and mass + additional units")) # unsupported
end
end
end

function kz{T0, T1, T2, T3, T4}(self::Atmosphere, df::T0, constant_kz::T1, chf::T2, kz_min::T3, latent_heat::T4)::String
return "Depricating this function. Please use ptk instead. It has identical functionality."
if !(isinstance(df, type_(nothing)))
self.chf = nothing
if df.loc[df["kz"] < kz_min].shape[0] > 0
df.loc[df["kz"] < kz_min] = kz_min
println(join(["Overwriting some Kz values to minimum value set by kz_min"], " "));
end
self.kz = array(np, df["kz"])
if length(self.kz) != length(self.pressure)
raise!(Exception("Kzz and pressure are not the same length")) # unsupported
end
else

if !(isinstance(constant_kz, type_(nothing)))
self.chf = nothing
self.kz = constant_kz
if self.kz < kz_min
self.kz = kz_min
println(join(["Overwriting kz constant value to minimum value set by kz_min"], " "));
end
else

if !(isinstance(chf, type_(nothing)))
function g_c_85{T0, T1, T2, T3, T4, T5, RT}(scale_h::T0, r_atmos::T1, chf::T2, rho_atmos::T3, c_p::T4, lapse_ratio::T5)::RT
if latent_heat
mixl = (max(np, 0.1, lapse_ratio)*scale_h)
else

mixl = scale_h
end
gc_kzz = ((((1.0/3.0)*scale_h)*pow((mixl/scale_h), (4.0/3.0)))*pow(((r_atmos*chf)/(rho_atmos*c_p)), (1.0/3.0)))
return (max(np, gc_kzz, kz_min), mixl)
end

self.kz = g_c_85
self.chf = chf
end
end
end
end

function compute{T0, T1, RT}(self::Atmosphere, directory::T0, as_dict::T1)::RT
run = compute(self)
return run
end

function calc_mie_db{T0, T1, T2, T3, T4, RT}(gas_name::T0, dir_refrind::T1, dir_out::T2, rmin::T3, nradii::T4)::RT
if isinstance(gas_name, str)
gas_name = [gas_name];
end
ngas = length(gas_name)
for i in (0:length(gas_name) - 1)
println(join([("Computing " + gas_name[i])], " "));
wave_in, nn, kk = get_refrind(gas_name[i], dir_refrind)
nwave = length(wave_in)
if i == 0
radius, rup, dr = get_r_grid()
qext_all = zeros(np, (nwave, nradii, ngas))
qscat_all = zeros(np, (nwave, nradii, ngas))
cos_qscat_all = zeros(np, (nwave, nradii, ngas))
end
qext_gas, qscat_gas, cos_qscat_gas = calc_new_mieff(wave_in, nn, kk, radius, rup, false)
qext_all[None], qscat_all[None], cos_qscat_all[None] = (qext_gas, qscat_gas, cos_qscat_gas)
wave = ([nwave] + radius.iter().map(|r| ([r] + list(wave_in))).collect::<Vec<_>>().iter().sum())
qscat = ([nradii] + qscat_gas.T.iter().map(|iscat| ([np.nan] + list(iscat))).collect::<Vec<_>>().iter().sum())
qext = ([np.nan] + qext_gas.T.iter().map(|iext| ([np.nan] + list(iext))).collect::<Vec<_>>().iter().sum())
cos_qscat = ([np.nan] + cos_qscat_gas.T.iter().map(|icos| ([np.nan] + list(icos))).collect::<Vec<_>>().iter().sum())
to_csv(pd.DataFrame(Dict("wave" => wave, "qscat" => qscat, "qext" => qext, "cos_qscat" => cos_qscat)), os.path.join(dir_out, (gas_name[i] + ".mieff")), " ", false, nothing);
end
return (qext_all, qscat_all, cos_qscat_all, radius, wave_in)
end

function get_mie{T0, T1, RT}(gas::T0, directory::T1)::RT
df = read_csv(pd, os.path.join(directory, (gas + ".mieff")), ["wave", "qscat", "qext", "cos_qscat"], true)
nwave = Int64(df.iloc[(0, 0)])
nradii = Int64(df.iloc[(0, 1)])
radii = df.loc[isnan(np, df["qscat"])]["wave"].values
df = dropna(df);
@assert(length(radii) == nradii)
@assert((nwave*nradii) == df.shape[0])
wave = reshape(None.values, (nradii, nwave)).T
qscat = reshape(None.values, (nradii, nwave)).T
qext = reshape(None.values, (nradii, nwave)).T
cos_qscat = reshape(None.values, (nradii, nwave)).T
return (qext, qscat, cos_qscat, nwave, radii, wave)
end

function get_refrind{T0, T1, RT}(igas::T0, directory::T1)::RT
filename = join(os.path, directory, (igas + ".refrind"))
idummy, wave_in, nn, kk = loadtxt(np, OpenOptions::new().read(true).open(filename).readlines(), true, [0, 1, 2, 3])
return (wave_in, nn, kk)
end

function get_r_grid_w_max{T0, T1, T2, RT}(r_min::T0, r_max::T1, n_radii::T2)::RT
radius = logspace(np, np.log10(r_min), np.log10(r_max), n_radii)
rat = (radius[1]/radius[0])
rup = (((2*rat)/(rat + 1))*radius)
dr = zeros(np, rup.shape)
dr[1..] = (rup[1..] - rup[..-1])
dr[0] = (pow(dr[1], 2)/dr[2])
return (radius, rup, dr)
end

function get_r_grid{T0, T1, RT}(r_min::T0, n_radii::T1)::RT
vrat = 2.2
pw = (1.0/3.0)
f1 = pow(((2.0*vrat)/(1.0 + vrat)), pw)
f2 = (pow((2.0/(1.0 + vrat)), pw)*pow(vrat, (pw - 1.0)))
radius = (r_min*pow(vrat, (linspace(np, 0, (n_radii - 1), n_radii)/3.0)))
rup = (f1*radius)
dr = (f2*radius)
return (radius, rup, dr)
end

function picaso_format{T0, T1, T2, RT}(opd::T0, w0::T1, g0::T2)::RT
df = DataFrame(pd, (0:(opd.shape[0]*opd.shape[1]) - 1).iter().map(|i| i).collect::<Vec<_>>(), ["lvl", "w", "opd", "w0", "g0"])
i = 0
LVL = []
WV, OPD, WW0, GG0 = ([], [], [], [])
for j in (0:opd.shape[0] - 1)
for w in (0:opd.shape[1] - 1)
LVL += [(j + 1)]
WV += [(w + 1)]
OPD += [opd[(j, w)]]
WW0 += [w0[(j, w)]]
GG0 += [g0[(j, w)]]
end
end
df.iloc[None] = LVL
df.iloc[None] = WV
df.iloc[None] = OPD
df.iloc[None] = WW0
df.iloc[None] = GG0
return df
end

function available{RT}()::RT
pvs = dir(pvaps).iter().cloned().filter(|&i| i != "np"&&"_" not in i).map(|i| i).collect::<Vec<_>>()
gas_p = dir(gas_properties).iter().cloned().filter(|&i| i != "np"&&"_" not in i).map(|i| i).collect::<Vec<_>>()
return list(np.intersect1d(gas_p, pvs))
end

function recommend_gas{T0, T1, T2, T3, T4, T5}(pressure::T0, temperature::T1, mh::T2, mmw::T3, plot::T4, legend::T5)::List
if plot
using bokeh::plotting: figure, show
using bokeh::models: Legend
using bokeh::palettes: magma
plot_kwargs["y_range"] = get(plot_kwargs, "y_range", [100.0, 0.001])
plot_kwargs["plot_height"] = get(plot_kwargs, "plot_height", 400)
plot_kwargs["plot_width"] = get(plot_kwargs, "plot_width", 600)
plot_kwargs["x_axis_label"] = get(plot_kwargs, "x_axis_label", "Temperature (K)")
plot_kwargs["y_axis_label"] = get(plot_kwargs, "y_axis_label", "Pressure (bars)")
plot_kwargs["y_axis_type"] = get(plot_kwargs, "y_axis_type", "log")
fig = figure(plot_kwargs)
end
all_gases = available()
cond_ts = []
recommend = []
line_widths = []
for gas_name in all_gases
cond_p, t = condensation_t(gas_name, mh, mmw)
cond_ts += [t]
interp_cond_t = interp(np, pressure, cond_p, t)
diff_curve = (interp_cond_t - temperature)
if (length(diff_curve[diff_curve > 0]) > 0 & length(diff_curve[diff_curve < 0]) > 0)
recommend += [gas_name]
line_widths += [5]
else

line_widths += [1]
end
end
if plot
legend_it = []
ngas = length(all_gases)
cols = magma(ngas)
if legend == "inside"
line(fig, temperature, pressure, "User", "black", 5, "dashed");
for i in (0:ngas - 1)
line(fig, cond_ts[i], cond_p, all_gases[i], cols[i], line_widths[i]);
end
else

f = line(fig, temperature, pressure, "black", 5, "dashed")
push!(legend_it, ("input profile", [f]));
for i in (0:ngas - 1)
f = line(fig, cond_ts[i], cond_p, cols[i], line_widths[i]);
push!(legend_it, (all_gases[i], [f]));
end
end
if legend == "outside"
legend = Legend(legend_it, (0, 0));
legend.click_policy = "mute"
add_layout(fig, legend, "right");
end
plot_format(fig);
show(fig);
end
return recommend
end

function condensation_t{T0, T1, T2, T3, RT}(gas_name::T0, mh::T1, mmw::T2, pressure::T3)::RT
if isinstance(pressure, (float, int))
pressure = [pressure];
end
temps = []
for p in pressure
temp = root_scalar(optimize, find_cond_t, [10, 10000], "brentq", (p, mh, mmw, gas_name))
temps += [temp.root]
end
return (array(np, pressure), array(np, temps))
end

function hot_jupiter{RT}()::RT
directory = join(os.path, os.path.dirname(__file__), "reference", "hj.pt")
df = read_csv(pd, directory, true, [1, 2, 3], ["pressure", "temperature", "kz"], 1)
df.loc[(df["pressure"] > 12.8, "temperature")] = linspace(np, 1822, 2100, df.loc[df["pressure"] > 12.8].shape[0])
return df
end

function brown_dwarf{RT}()::RT
directory = join(os.path, os.path.dirname(__file__), "reference", "t1000g100nc_m0.0.dat")
df = read_csv(pd, directory, 1, true, nothing, [1, 2, 3], ["pressure", "temperature", "chf"])
return df
end

