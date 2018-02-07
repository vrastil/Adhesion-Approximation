
constexpr FTYPE PI = M_PI;
constexpr FTYPE MPL = 2.435E18; // reduced Planck mass, [GeV/c^2]
constexpr FTYPE fm_to_Mpc = 3.2407792896664E-38; // 1 fm = ? Mpc
constexpr FTYPE hbarc = 197.327053; // reduced Planck constant times speed of light, [MeV fm]
constexpr FTYPE hbarc_cosmo = hbarc*1E-9*fm_to_Mpc; // [GeV Mpc]
constexpr FTYPE G_N = hbarc_cosmo*6.70711*1E-39*hbarc_cosmo; // gravitational constant, [GeV Mpc / (GeV/c^2)^2]
constexpr FTYPE c_kms = 299792.458; // speed of light [km / s]