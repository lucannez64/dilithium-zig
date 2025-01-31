const poly = @import("poly.zig");
const params = @import("params.zig");
pub const polyvecl = struct {
    vec: [params.L]poly.poly,
};

pub const polyveck = struct {
    vec: [params.K]poly.poly,
};

pub fn polyvec_matrix_expand(mat: *[params.K]polyvecl, rho: [params.SEEDBYTES]u8) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        var j: usize = 0;
        while (j < params.L) : (j += 1) {
            poly.poly_uniform(&mat[i].vec[j], rho, (i << 8) + j);
        }
    }
}

pub fn polyvec_matrix_pointwise_montgomery(t: *polyveck, mat: [params.K]polyvecl, v: *polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        polyvecl_pointwise_acc_montgomery(&t.vec[i], &mat[i].vec, v);
    }
}

pub fn polyvecl_uniform_eta(v: *polyvecl, seed: [params.CRHBYTES]u8, nonce: u16) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_uniform_eta(&v.vec[i], seed, nonce);
        nonce += 1;
    }
}

pub fn polyvecl_uniform_gamma1(v: *polyvecl, seed: [params.CRHBYTES]u8, nonce: u16) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_uniform_gamma1(&v.vec[i], seed, params.L * nonce + i);
        nonce += 1;
    }
}

pub fn polyvecl_reduce(v: *polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_reduce(&v.vec[i]);
    }
}

pub fn polyvecl_add(w: *polyvecl, u: *const polyvecl, v: *const polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_add(&w.vec[i], &u.vec[i], &v.vec[i]);
    }
}

pub fn polyvecl_ntt(v: *polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_ntt(&v.vec[i]);
    }
}

pub fn polyvecl_invntt_tomont(v: *polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_invntt_tomont(&v.vec[i]);
    }
}

pub fn polyvecl_pointwise_poly_montgomery(r: *polyvecl, a: *const poly.poly, v: *const polyvecl) void {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        poly.poly_pointwise_montgomery(&r.vec[i], a, &v.vec[i]);
    }
}

pub fn polyvecl_pointwise_acc_montgomery(r: *poly, a: *const polyvecl, v: *const polyvecl) void {
    var i: usize = 0;
    var t: poly.poly = undefined;
    poly.poly_pointwise_montgomery(r, &a.vec[0], &v.vec[0]);
    while (i < params.L) : (i += 1) {
        poly.poly_pointwise_acc_montgomery(&t, &a.vec[i], &v.vec[i]);
        poly.poly_add(r, r, &t);
    }
}

pub fn polyvecl_chknorm(v: *const polyvecl, bound: i32) bool {
    var i: usize = 0;
    while (i < params.L) : (i += 1) {
        if (poly.poly_chknorm(&v.vec[i], bound)) {
            return true;
        }
    }
    return false;
}

pub fn polyveck_uniform_eta(v: *polyveck, seed: [params.CRHBYTES]u8, nonce: u16) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_uniform_eta(&v.vec[i], seed, nonce);
        nonce += 1;
    }
}

pub fn polyveck_reduce(v: *polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_reduce(&v.vec[i]);
    }
}

pub fn polyveck_caddq(v: *polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_caddq(&v.vec[i]);
    }
}

pub fn polyveck_add(w: *polyveck, u: *const polyveck, v: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_add(&w.vec[i], &u.vec[i], &v.vec[i]);
    }
}

pub fn polyveck_sub(w: *polyveck, u: *const polyveck, v: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_sub(&w.vec[i], &u.vec[i], &v.vec[i]);
    }
}

pub fn polyveck_shiftl(v: *polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_shiftl(&v.vec[i]);
    }
}

pub fn polyveck_ntt(v: *polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_ntt(&v.vec[i]);
    }
}

pub fn polyveck_invntt_tomont(v: *polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_invntt_tomont(&v.vec[i]);
    }
}

pub fn polyveck_pointwise_poly_montgomery(r: *polyveck, a: *const poly.poly, v: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_pointwise_montgomery(&r.vec[i], a, &v.vec[i]);
    }
}

pub fn polyveck_chknorm(v: *const polyveck, bound: i32) bool {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        if (poly.poly_chknorm(&v.vec[i], bound)) {
            return true;
        }
    }
    return false;
}

pub fn polyveck_power2round(v1: *polyveck, v0: *polyveck, v: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_power2round(&v1.vec[i], &v0.vec[i], &v.vec[i]);
    }
}

pub fn polyveck_decompose(v1: *polyveck, v0: *polyveck, v: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_decompose(&v1.vec[i], &v0.vec[i], &v.vec[i]);
    }
}

pub fn polyveck_make_hint(h: *polyveck, v0: *const polyveck, v1: *const polyveck) usize {
    var i: usize = 0;
    var s: usize = 0;
    while (i < params.K) : (i += 1) {
        s += poly.poly_make_hint(&h.vec[i], &v0.vec[i], &v1.vec[i]);
    }
    return s;
}

pub fn polyveck_use_hint(w: *polyveck, u: *const polyveck, h: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_use_hint(&w.vec[i], &u.vec[i], &h.vec[i]);
    }
}

pub fn polyveck_pack_w1(r: [params.K * params.POLYW1_PACKEDBYTES]u8, w1: *const polyveck) void {
    var i: usize = 0;
    while (i < params.K) : (i += 1) {
        poly.poly_pack_w1(&r[i * params.POLYW1_PACKEDBYTES], &w1.vec[i]);
    }
}
