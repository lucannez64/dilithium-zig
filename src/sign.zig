const params = @import("params.zig");
const std = @import("std");
const poly = @import("poly.zig");
const polyvec = @import("polyvec.zig");
const packing = @import("packing.zig");
const symmetric = @import("symmetric.zig");
const SEEDBYTES = params.SEEDBYTES;
const TRBYTES = params.TRBYTES;
const CRHBYTES = params.CRHBYTES;
fn randombytes(buf: []u8) void {
    std.crypto.random.bytes(buf);
}

pub fn crypto_sign_keypair(pk: []u8, sk: []u8) void {
    var buffer: [SEEDBYTES + TRBYTES + SEEDBYTES + CRHBYTES]u8 = undefined;
    var seedbuf: [2 * params.SEEDBYTES + params.CRHBYTES]u8 = buffer[0 .. 2 * params.SEEDBYTES + params.CRHBYTES];
    var tr: [params.TRBYTES]u8 = undefined;
    var rho: []u8 = undefined;
    var rhoprime: []u8 = undefined;
    var key: []u8 = undefined;
    var mat: [params.K]polyvec.polyvecl = undefined;
    var s1: polyvec.polyvecl = undefined;
    var s1hat: polyvec.polyvecl = undefined;
    var s2: polyvec.polyveck = undefined;
    var t1: polyvec.polyveck = undefined;
    var t0: polyvec.polyveck = undefined;

    randombytes(&seedbuf);
    seedbuf[params.SEEDBYTES] = params.K;
    seedbuf[params.SEEDBYTES + 1] = params.L;
    const seed = seedbuf[0 .. params.SEEDBYTES + 2];
    std.crypto.hash.sha3.Shake256.hash(seed, &seedbuf, .{});
    rho = buffer[0..params.SEEDBYTES];
    rhoprime = buffer[params.SEEDBYTES .. params.SEEDBYTES + params.CRHBYTES];
    key = buffer[params.SEEDBYTES + params.CRHBYTES .. params.SEEDBYTES + params.CRHBYTES + params.SEEDBYTES];
    polyvec.polyvec_matrix_expand(&mat, &rho);
    polyvec.polyvecl_uniform_eta(&s1, rhoprime, 0);
    polyvec.polyveck_uniform_eta(&s2, rhoprime, params.L);
    s1hat = s1;
    polyvec.polyvecl_ntt(&s1hat);

    polyvec.polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
    polyvec.polyveck_reduce(&t1);
    polyvec.polyveck_invntt_tomont(&t1);

    polyvec.polyveck_add(&t1, &t1, &s2);
    polyvec.polyveck_caddq(&t1);
    polyvec.polyveck_power2round(&t1, &t0, &t1);
    packing.pack(pk, rho, &t1);

    polyvec.polyveck_caddq(&t1);
    polyvec.polyveck_power2round(&t1, &t0, &t1);
    packing.pack_pk(pk[0..params.CRYPTO_PUBLICKEYBYTES].*, rho[0..params.SEEDBYTES].*, &t1, &t1);
    std.crypto.hash.sha3.Shake256.hash(pk, &tr, .{});
    packing.pack_sk(sk[0..params.CRYPTO_SECRETKEYBYTES].*, rho[0..params.SEEDBYTES].*, tr[0..params.TRBYTES].*, key[0..params.SEEDBYTES].*, &t0, &s1, &s2);
}

pub fn crypto_sign_signature_internal(
    sig: *[]u8,
    siglen: *usize,
    m: []const u8,
    pre: []u8,
    rnd: []const u8,
    sk: []const u8,
) isize {
    var buffer: [SEEDBYTES + TRBYTES + SEEDBYTES + CRHBYTES + CRHBYTES]u8 = undefined;
    var rho: []u8 = undefined;
    var tr: []u8 = undefined;
    var key: []u8 = undefined;
    var mu: []u8 = undefined;
    var rhoprime: [params.CRHBYTES]u8 = undefined;
    var nonce: u16 = 0;
    var mat: [params.K]polyvec.polyvecl = undefined;
    var s1: polyvec.polyvecl = undefined;
    var y: polyvec.polyvecl = undefined;
    var z: polyvec.polyvecl = undefined;
    var t0: polyvec.polyveck = undefined;
    var s2: polyvec.polyveck = undefined;
    var w1: polyvec.polyveck = undefined;
    var w0: polyvec.polyveck = undefined;
    var h: polyvec.polyveck = undefined;
    var cp: poly.poly = undefined;
    var ctx = symmetric.init256();
    var state = ctx.state.st;
    rho = &buffer[0..params.SEEDBYTES].*;
    tr = buffer[params.SEEDBYTES .. params.SEEDBYTES + params.TRBYTES];
    key = buffer[params.SEEDBYTES + params.TRBYTES .. params.SEEDBYTES + params.TRBYTES + params.SEEDBYTES];
    mu = buffer[params.SEEDBYTES + params.TRBYTES + params.SEEDBYTES .. params.SEEDBYTES + params.TRBYTES + params.SEEDBYTES + params.CRHBYTES];
    rhoprime = buffer[params.SEEDBYTES + params.TRBYTES + params.SEEDBYTES + params.CRHBYTES .. params.SEEDBYTES + params.TRBYTES + params.SEEDBYTES + params.CRHBYTES + params.CRHBYTES].*;

    packing.unpack_sk(rho[0..params.SEEDBYTES], tr[0..params.TRBYTES], key[0..params.SEEDBYTES], &t0, &s1, &s2, sk[0..params.CRYPTO_SECRETKEYBYTES]);

    // Compute mu = CRH(tr, pre, msg)
    state.absorb(tr);
    state.absorb(pre);
    state.absorb(m);
    ctx.final();
    state.squeeze(mu);

    // Compute rhoprime = CRH(key, rnd, mu)
    ctx = symmetric.init256();
    state = ctx.state.st;
    state.absorb(key);
    state.absorb(rnd);
    state.absorb(mu);
    ctx.final();
    state.squeeze(&rhoprime);

    // Expand matrix and transform vectors
    polyvec.polyvec_matrix_expand(mat, &rho);
    polyvec.polyvecl_ntt(&s1);
    polyvec.polyveck_ntt(&s2);
    polyvec.polyveck_ntt(&t0);

    rej: while (true) {
        // Sample intermediate vector y
        polyvec.polyvecl_uniform_gamma1(&y, rhoprime, nonce);
        nonce += 1;

        // Matrix-vector multiplication
        z = y;
        polyvec.polyvecl_ntt(&z);
        polyvec.polyvec_matrix_pointwise_montgomery(&w1, &mat, &z);
        polyvec.polyveck_reduce(&w1);
        polyvec.polyveck_invntt_tomont(&w1);

        // Decompose w and call the random oracle
        polyvec.polyveck_caddq(&w1);
        polyvec.polyveck_decompose(&w1, &w0, &w1);
        polyvec.polyveck_pack_w1(sig, &w1);

        ctx = symmetric.init256();
        state = ctx.state.st;
        state.absorb(mu);
        state.absorb(sig[0 .. params.K * params.POLYW1_PACKEDBYTES]);
        ctx.final();
        state.squeeze(&sig[0..params.CTILDEBYTES]);

        poly.poly_challenge(&cp, sig);
        poly.poly_ntt(&cp);

        // Compute z, reject if it reveals secret
        polyvec.polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
        polyvec.polyvecl_invntt_tomont(&z);
        polyvec.polyvecl_add(&z, &z, &y);
        polyvec.polyvecl_reduce(&z);
        if (polyvec.polyvecl_chknorm(&z, params.GAMMA1 - params.BETA)) {
            continue :rej;
        }

        // Check that subtracting cs2 does not change high bits of w and low bits
        // do not reveal secret information
        polyvec.polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
        polyvec.polyveck_invntt_tomont(&h);
        polyvec.polyveck_sub(&w0, &w0, &h);
        polyvec.polyveck_reduce(&w0);
        if (polyvec.polyveck_chknorm(&w0, params.GAMMA2 - params.BETA)) {
            continue :rej;
        }

        // Compute hints for w1
        polyvec.polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
        polyvec.polyveck_invntt_tomont(&h);
        polyvec.polyveck_reduce(&h);
        if (polyvec.polyveck_chknorm(&h, params.GAMMA2)) {
            continue :rej;
        }

        polyvec.polyveck_add(&w0, &w0, &h);
        const n = polyvec.polyveck_make_hint(&h, &w0, &w1);
        if (n > params.OMEGA) {
            continue :rej;
        }

        // Write signature
        packing.pack_sig(sig, sig, &z, &h);
        siglen.* = params.CRYPTO_BYTES;
        return 0;
    }
}

pub const SignatureError = error{
    ContextTooLong,
    BadSignature,
};

pub fn crypto_sign_signature(sig: *[]u8, m: []const u8, ctx: []const u8, sk: []const u8) SignatureError!isize {
    var i: usize = 0;
    var pre: [256]u8 = undefined;
    var rnd: [params.RNDBYTES]u8 = undefined;
    if (ctx.len > 256) {
        return SignatureError.ContextTooLong;
    }
    pre[0] = 0;
    pre[1] = @truncate(ctx.len);
    while (i < ctx.len) : (i += 1) {
        pre[i + 2] = ctx[i];
    }
    randombytes(&rnd);
    return crypto_sign_signature_internal(sig, &sig.len, m, pre[0 .. 2 + ctx.len], &rnd, sk);
}

pub fn crypto_sign(sm: []u8, m: []u8, ctx: []const u8, sk: []const u8) SignatureError!isize {
    var i: usize = 0;
    while (i < m.len) : (i += 1) {
        sm[params.CRYPTO_BYTES + m.len - 1 - i] = m[m.len - 1 - i];
    }
    const ret = crypto_sign_signature(sm[0..], sm[params.CRYPTO_BYTES..], ctx, sk);
    return ret;
}

pub fn crypto_sign_verify_internal(sig: []const u8, m: []const u8, pre: []u8, pk: []const u8) bool {
    var i: usize = 0;
    var buf: [params.K * params.POLYW1_PACKEDBYTES]u8 = undefined;
    var rho: [params.SEEDBYTES]u8 = undefined;
    var mu: [params.CRHBYTES]u8 = undefined;
    var c: [params.CTILDEBYTES]u8 = undefined;
    var c2: [params.CTILDEBYTES]u8 = undefined;
    var cp: poly.poly = undefined;
    var mat: [params.K]polyvec.polyvecl = undefined;
    var z: polyvec.polyvecl = undefined;
    var t1: polyvec.polyveck = undefined;
    var w1: polyvec.polyveck = undefined;
    var h: polyvec.polyveck = undefined;
    var ctx: symmetric.Shake256_state = symmetric.init256();
    var state: std.crypto.hash.sha3.Shake256 = ctx.state;

    if (sig.len != params.CRYPTO_BYTES) {
        return false;
    }

    packing.unpack_pk(&rho, &t1, pk);
    if (packing.unpack_sig(c, &z, &h, sig)) {
        return false;
    }
    if (polyvec.polyvecl_chknorm(&z, params.GAMMA1 - params.BETA)) {
        return false;
    }
    state.hash(pk, mu[0..params.TRBYTES], .{});
    ctx = symmetric.init256();
    state = ctx.state.st;
    state.absorb(mu[0..params.TRBYTES]);
    state.absorb(pre);
    state.absorb(m);
    ctx.final();
    state.squeeze(&mu);
    poly.poly_challenge(&cp, &c);
    polyvec.polyvec_matrix_expand(&mat, &rho);
    polyvec.polyvecl_ntt(&z);
    polyvec.polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
    poly.poly_ntt(&cp);
    polyvec.polyveck_shiftl(&t1);
    polyvec.polyveck_ntt(&t1);
    polyvec.polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);
    polyvec.polyveck_sub(&w1, &w1, &t1);
    polyvec.polyveck_reduce(&w1);
    polyvec.polyveck_invntt_tomont(&w1);
    polyvec.polyveck_caddq(&w1);
    polyvec.polyveck_use_hint(&w1, &w1, &h);
    polyvec.polyveck_pack_w1(&buf, &w1);
    ctx = symmetric.init256();
    state = ctx.state.st;
    state.absorb(mu);
    state.absorb(buf);
    ctx.final();
    state.squeeze(&c2);
    i = 0;
    while (i < params.CTILDEBYTES) : (i += 1) {
        if (c[i] != c2[i]) {
            return false;
        }
    }
    return true;
}

pub fn crypto_sign_verify(sig: []const u8, m: []const u8, ctx: []const u8, pk: []const u8) SignatureError!bool {
    var i: usize = 0;
    var pre: [257]u8 = undefined;
    if (ctx.len > 255) {
        return SignatureError.ContextTooLong;
    }
    pre[0] = 0;
    pre[1] = @truncate(ctx.len);
    while (i < ctx.len) : (i += 1) {
        pre[i + 2] = ctx[i];
    }
    return crypto_sign_verify_internal(sig, m, pre[0 .. ctx.len + 2], pk);
}

pub fn crypto_sign_open(m: []u8, sm: []const u8, ctx: []const u8, pk: []const u8) SignatureError!void {
    var i: usize = 0;
    if (sm.len < params.CRYPTO_BYTES) {
        badsig(&m, sm.len);
        return SignatureError.BadSignature;
    }

    const mlen = sm.len - params.CRYPTO_BYTES;
    if (crypto_sign_verify(sm[0..params.CRYPTO_BYTES], sm[params.CRYPTO_BYTES .. params.CRYPTO_BYTES + mlen], ctx, pk)) {
        badsig(&m, sm.len);
        return SignatureError.BadSignature;
    } else {
        while (i < mlen) : (i += 1) {
            m[i] = sm[params.CRYPTO_BYTES + i];
        }
        return;
    }
}

fn badsig(m: []u8, smlen: usize) void {
    var i: usize = 0;
    while (i < smlen) : (i += 1) {
        m[i] = 0;
        m.len = 0;
    }
}
