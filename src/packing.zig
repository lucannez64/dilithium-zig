const polyvec = @import("polyvec.zig");
const poly = @import("poly.zig");
const params = @import("params.zig");

pub fn pack_pk(pk: [params.CRYPTO_PUBLICKEYBYTES]u8, rho: [params.SEEDBYTES]u8, t1: *const polyvec.polyveck) void {
    var i: usize = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        pk[i] = rho[i];
    }
    pk += params.SEEDBYTES;

    i = 0;
    while (i < params.POLYVECK_N) : (i += 1) {
        polyvec.pack(pk + i * params.POLYT1_PACKEDBYTES, &t1.vec[i]);
    }
}

pub fn unpack_pk(rho: [params.SEEDBYTES]u8, t1: *polyvec.polyveck, pk: [params.CRYPTO_PUBLICKEYBYTES]u8) void {
    var i: usize = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        rho[i] = pk[i];
    }
    pk += params.SEEDBYTES;

    i = 0;
    while (i < params.POLYVECK_N) : (i += 1) {
        polyvec.unpack(&t1.vec[i], pk + i * params.POLYT1_PACKEDBYTES);
    }
}

pub fn pack_sk(sk: *[params.CRYPTO_SECRETKEYBYTES]u8, rho: [params.SEEDBYTES]u8, tr: [params.TRBYTES]u8, key: [params.SEEDBYTES]u8, t0: *const polyvec.polyveck, s1: *const polyvec.polyvecl, s2: *const polyvec.polyveck) void {
    var i: usize = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        sk[i] = rho[i];
    }

    i = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        sk[i + params.SEEDBYTES] = key[i];
    }

    i = 0;
    while (i < params.TRBYTES) : (i += 1) {
        sk[i + 2 * params.SEEDBYTES] = tr[i];
    }

    i = 0;
    while (i < params.L) : (i += 1) {
        poly.polyeta_pack(sk[(i + 2 * params.SEEDBYTES + params.TRBYTES) * params.POLYETA_PACKEDBYTES ..].*, &s1.vec[i]);
    }

    i = 0;
    while (i < params.K) : (i += 1) {
        polyvec.pack(sk[(i + 2 * params.SEEDBYTES + params.TRBYTES + params.L * params.POLYETA_PACKEDBYTES) * params.POLYT0_PACKEDBYTES ..].*, &s2.vec[i]);
    }

    i = 0;
    while (i < params.K) : (i += 1) {
        polyvec.pack(sk[(i + 2 * params.SEEDBYTES + params.TRBYTES + params.L * params.POLYETA_PACKEDBYTES + params.K * params.POLYT0_PACKEDBYTES) * params.POLYT0_PACKEDBYTES ..].*, &t0.vec[i]);
    }
}

pub fn unpack_sk(rho: *[params.SEEDBYTES]u8, tr: *[params.TRBYTES]u8, key: *[params.SEEDBYTES]u8, t0: *polyvec.polyveck, s1: *polyvec.polyvecl, s2: *polyvec.polyveck, sk: *const [params.CRYPTO_SECRETKEYBYTES]u8) void {
    var i: usize = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        rho[i] = sk[i];
    }

    var j: usize = i;
    i = 0;
    while (i < params.SEEDBYTES) : (i += 1) {
        key[i] = sk[j];
        j += 1;
    }

    i = 0;
    while (i < params.TRBYTES) : (i += 1) {
        tr[i] = sk[j];
        j += 1;
    }

    i = 0;
    while (i < params.L) : (i += 1) {
        poly.polyeta_unpack(&s1.vec[i], sk[j * params.POLYT0_PACKEDBYTES ..]);
        j += 1;
    }

    i = 0;
    while (i < params.K) : (i += 1) {
        poly.polyeta_unpack(&s2.vec[i], sk[j * params.POLYT0_PACKEDBYTES ..]);
        j += 1;
    }

    i = 0;
    while (i < params.K) : (i += 1) {
        poly.polyt0_unpack(&t0.vec[i], sk[j * params.POLYT0_PACKEDBYTES ..]);
        j += 1;
    }
}

pub fn pack_sig(sig: [params.CRYPTO_BYTES]u8, c: [params.CTILDEBYTES]u8, z: *const polyvec.polyvecl, h: *const polyvec.polyveck) void {
    var i: usize = 0;
    while (i < params.CTILDEBYTES) : (i += 1) {
        sig[i] = c[i];
    }

    var eee: usize = i;
    i = 0;
    while (i < params.L) : (i += 1) {
        poly.polyeta_pack(sig[eee * params.POLYETA_PACKEDBYTES ..].*, &z.vec[i]);
        eee += 1;
    }

    i = 0;
    while (i < params.OMEGA + params.K) : (i += 1) {
        sig[eee] = 0;
        eee += 1;
    }
    var k: usize = 0;
    i = 0;
    while (i < params.K) : (i += 1) {
        var j = 0;
        while (j < params.N) : (j += 1) {
            if (h.vec[i].coeffs[j] != 0) {
                sig[k] = j;
                k += 1;
            }
        }
        sig[params.OMEGA + i] = k;
    }
}

pub fn unpack_sig(c: [params.CTILDEBYTES]u8, z: *polyvec.polyvecl, h: *polyvec.polyveck, sig: [params.CRYPTO_BYTES]u8) bool {
    var i: usize = 0;
    while (i < params.CTILDEBYTES) : (i += 1) {
        c[i] = sig[i];
    }
    var eee: usize = i;
    i = 0;
    while (i < params.L) : (i += 1) {
        poly.polyeta_unpack(&z.vec[i], sig[eee * params.POLYETA_PACKEDBYTES ..]);
        eee += 1;
    }
    sig = sig[params.L * params.POLYETA_PACKEDBYTES ..].*;

    var k: usize = 0;
    i = 0;

    var j: usize = 0;
    while (i < params.K) : (i += 1) {
        j = 0;
        while (j < params.N) : (j += 1) {
            h.vec[i].coeffs[j] = 0;
        }
        if ((sig[params.OMEGA + i] < k) or (sig[params.OMEGA + i] > params.OMEGA)) {
            return true;
        }
        j = k;
        while (j < sig[params.OMEGA + i]) : (j += 1) {
            if ((j > k) and (sig[j] <= sig[j - 1])) {
                return true;
            }
            h.vec[i].coeffs[sig[j]] = 1;
        }

        k = sig[params.OMEGA + i];
    }

    j = k;

    while (j < params.OMEGA) : (j += 1) {
        if (sig[j]) {
            return true;
        }
    }
    return 0;
}
