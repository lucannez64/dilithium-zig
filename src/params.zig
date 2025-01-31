pub const config = @import("config.zig");
pub const DILITHIUM_MODE = config.DILITHIUM_MODE;
pub const RNDBYTES = 32;
pub const SEEDBYTES = 32;
pub const CRHBYTES = 64;
pub const TRBYTES = 64;
pub const N = 256;
pub const Q = 8380417;
pub const D = 13;
pub const ROOT_OF_UNITY = 1753;

pub const K: usize = if (DILITHIUM_MODE == 2) 4 else if (DILITHIUM_MODE == 3) 6 else 8;
pub const L: usize = if (DILITHIUM_MODE == 2) 4 else if (DILITHIUM_MODE == 3) 5 else 7;
pub const ETA: usize = if (DILITHIUM_MODE == 2) 2 else if (DILITHIUM_MODE == 3) 4 else 2;
pub const TAU: usize = if (DILITHIUM_MODE == 2) 39 else if (DILITHIUM_MODE == 3) 49 else 60;
pub const BETA: usize = if (DILITHIUM_MODE == 2) 78 else if (DILITHIUM_MODE == 3) 196 else 120;
pub const GAMMA1: usize = if (DILITHIUM_MODE == 2) (1 << 17) else (1 << 19);
pub const GAMMA2: usize = if (DILITHIUM_MODE == 2) ((Q - 1) / 88) else ((Q - 1) / 32);
pub const OMEGA: usize = if (DILITHIUM_MODE == 2) 80 else if (DILITHIUM_MODE == 3) 55 else 75;
pub const CTILDEBYTES: usize = if (DILITHIUM_MODE == 2) 32 else if (DILITHIUM_MODE == 3) 48 else 64;

pub const POLYT1_PACKEDBYTES = 320;
pub const POLYT0_PACKEDBYTES = 416;
pub const POLYVECH_PACKEDBYTES = OMEGA + K;

pub const POLYZ_PACKEDBYTES = if (GAMMA1 == (1 << 17)) 576 else if (GAMMA1 == (1 << 19)) 640 else unreachable;

pub const POLYW1_PACKEDBYTES = if (GAMMA2 == ((Q - 1) / 88)) 192 else if (GAMMA2 == ((Q - 1) / 32)) 128 else unreachable;

pub const POLYETA_PACKEDBYTES = if (ETA == 2) 96 else if (ETA == 4) 128 else unreachable;

pub const CRYPTO_PUBLICKEYBYTES = SEEDBYTES + K * POLYT1_PACKEDBYTES;
pub const CRYPTO_SECRETKEYBYTES = 2 * SEEDBYTES +
    TRBYTES +
    L * POLYETA_PACKEDBYTES +
    K * POLYETA_PACKEDBYTES +
    K * POLYT0_PACKEDBYTES;
pub const CRYPTO_BYTES = CTILDEBYTES + L * POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES;
