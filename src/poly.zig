const params = @import("params.zig");
const ntt = @import("ntt.zig");
const std = @import("std");
const reduce = @import("reduce.zig");
const rounding = @import("rounding.zig");
const symmetric = @import("symmetric.zig");

pub const poly = struct {
    coeffs: [params.N]i32,
};

pub fn poly_reduce(p: *poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        p.coeffs[i] = reduce.reduce32(p.coeffs[i]);
    }
}

pub fn poly_caddq(p: *poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        p.coeffs[i] = reduce.caddq(p.coeffs[i]);
    }
}

pub fn poly_add(c: *poly, a: *const poly, b: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        c.coeffs[i] = a.coeffs[i] + b.coeffs[i];
    }
}

pub fn poly_sub(c: *poly, a: *const poly, b: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        c.coeffs[i] = a.coeffs[i] - b.coeffs[i];
    }
}

pub fn poly_shiftl(a: *poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        a.coeffs[i] <<= params.D;
    }
}

pub fn poly_ntt(a: *poly) void {
    ntt.ntt(a.coeffs);
}

pub fn poly_invntt_tomont(a: *poly) void {
    ntt.invntt_tomont(a.coeffs);
}

pub fn poly_pointwise_montgomery(c: *poly, a: *const poly, b: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        c.coeffs[i] = reduce.montgomery_reduce(@as(i64, a.coeffs[i]) * b.coeffs[i]);
    }
}

pub fn poly_power2round(a1: *poly, a0: *poly, a: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        a1.coeffs[i] = rounding.power2round(&a0.coeffs[i], a.coeffs[i]);
    }
}

pub fn poly_decompose(a1: *poly, a0: *poly, a: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        a1.coeffs[i] = symmetric.decompose(&a0.coeffs[i], a.coeffs[i]);
    }
}

pub fn poly_make_hint(h: *poly, a0: *poly, a1: *const poly) usize {
    var i: usize = 0;
    var s: usize = 0;
    while (i < params.N) : (i += 1) {
        h.coeffs[i] = symmetric.make_hint(a0.coeffs[i], a1.coeffs[i]);
        s += h.coeffs[i];
    }

    return s;
}

pub fn poly_use_hint(b: *poly, a: *const poly, h: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        b.coeffs[i] = symmetric.use_hint(a.coeffs[i], h.coeffs[i]);
    }
}

pub fn poly_chknorm(a: *const poly, B: i32) bool {
    var i: usize = 0;
    var t: i32 = 0;
    if (B > (params.Q - 1) / 8) {
        return true;
    }
    while (i < params.N) : (i += 1) {
        t = a.coeffs[i] >> 31;
        t = a.coeffs[i] - (t & 2 * a.coeffs[i]);
        if (t >= B) {
            return true;
        }
    }
    return false;
}

pub fn rej_uniform(a: *i32, len: usize, buf: *const u8, buflen: usize) usize {
    var ctr: usize = 0;
    var pos: usize = 0;
    var t: u32 = 0;
    while ((ctr < len) and (pos + 3 <= buflen)) {
        t = buf[pos];
        pos += 1;
        t |= @as(u32, buf[pos]) << 8;
        pos += 1;
        t |= @as(u32, buf[pos]) << 16;
        pos += 1;
        t &= "\x7FFFFF";
        if (t < params.Q) {
            a[ctr] = t;
            ctr += 1;
        }
    }

    return ctr;
}

const POLY_UNIFORM_NBLOCKS = ((768 + 168 - 1) / 168);
const STREAM128_BLOCKBYTES = 168;

pub fn poly_uniform(a: *poly, seed: [params.SEEDBYTES]u8, nonce: u16) void {
    var i = 0;
    var state = symmetric.init128();
    state.state.update(seed);

    var nonce_bytes: [2]u8 = undefined;
    std.mem.writeIntLittle(u16, &nonce_bytes, nonce);
    state.state.update(nonce_bytes);

    state.state.finalize();

    const total_block_size = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2;
    var buf: [total_block_size]u8 = undefined;

    state.state.squeeze(buf[0..(POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES)]);

    var ctr: usize = rej_uniform(a.coeffs[0..], params.N, &buf, POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES);

    while (ctr < params.N) {
        const off = (POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES) % 3;
        while (i < off) : (i += 1) {
            buf[i] = buf[POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES - off + i];
        }

        state.state.squeeze(buf[off..]);
        const buflen = STREAM128_BLOCKBYTES + off;
        _ = state.state.squeeze(buf[off..][0 .. buf.len - off]);
        ctr += rej_uniform(a.coeffs[ctr..], params.N - ctr, &buf, buflen);
    }
}

pub fn rej_eta(a: *i32, len: usize, buf: *const u8, buflen: usize) usize {
    var ctr: usize = 0;
    var pos: usize = 0;
    var t0: u32 = 0;
    var t1: u32 = 0;
    while ((ctr < len) and (pos < buflen)) {
        t0 = buf[pos] & "\x0F";
        t1 = buf[pos] >> 4;
        pos += 1;
    }

    if (params.ETA == 2) {
        if (t0 < 15) {
            t0 = t0 - (205 * t0 >> 10) * 5;
            a[ctr] = 2 - t0;
            ctr += 1;
        }
        if ((t1 < 15) and (ctr < len)) {
            t1 = t1 - (205 * t1 >> 10) * 5;
            a[ctr] = 2 - t1;
            ctr += 1;
        }
    } else if (params.ETA == 4) {
        if (t0 < 9) {
            a[ctr] = 4 - t0;
            ctr += 1;
        }
        if ((t1 < 9) and (ctr < len)) {
            a[ctr] = 4 - t1;
            ctr += 1;
        }
    }

    return ctr;
}

const STREAM256_BLOCKBYTES = 136;
const SHAKE256_RATE = 136;

const POLY_UNIFORM_ETA_NBLOCKS_2 = (227 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;
const POLY_UNIFORM_ETA_NBLOCKS_4 = (136 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;
const POLY_UNIFORM_ETA_NBLOCKS = if (params.ETA == 4) POLY_UNIFORM_ETA_NBLOCKS_4 else POLY_UNIFORM_ETA_NBLOCKS_2;

pub fn poly_uniform_eta(a: *poly, seed: [params.CRHBYTES]u8, nonce: u16) void {
    // Initialize SHAKE256 sponge (equivalent to stream256)
    var state = symmetric.init256();
    state.state.update(seed);

    // Process nonce into the state
    var nonce_bytes: [2]u8 = undefined;
    std.mem.writeInt(u16, &nonce_bytes, nonce, .little);
    state.state.update(nonce_bytes);
    state.state.finalize(); // Start squeezing output

    // Allocate buffer (same as C version)
    const buf_size = POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES;
    var buf: [buf_size]u8 = undefined;

    // Generate initial blocks
    state.state.squeeze(buf[0..buf_size]);

    var coefficients_filled: usize = 0;

    // Fill coefficients using initial buffer
    coefficients_filled = rej_eta(a.coeffs[0..], coefficients_filled, buf[0..], buf_size);

    // Continue generating blocks until all coefficients are filled
    while (coefficients_filled < params.N) {
        // Squeeze new block into the buffer (overwriting first block)
        _ = state.state.squeeze(buf[0..STREAM256_BLOCKBYTES]);

        // Fill remaining coefficients
        const added = rej_eta(a.coeffs[coefficients_filled..], params.N - coefficients_filled, buf[0..STREAM256_BLOCKBYTES], STREAM256_BLOCKBYTES);
        coefficients_filled += added;
    }
}

const POLY_UNIFORM_GAMMA1_NBLOCKS = ((params.POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES);

pub fn poly_uniform_gamma1(a: *poly, seed: [params.CRHBYTES]u8, nonce: u16) void {
    var state = symmetric.init256();
    var ctx = &state.state;
    ctx.update(seed);

    var nonce_bytes: [2]u8 = undefined;
    std.mem.writeIntLittle(u16, &nonce_bytes, nonce);
    ctx.update(nonce_bytes);

    var result: [POLY_UNIFORM_GAMMA1_NBLOCKS * 32]u8 = undefined;
    ctx.finalHash(result[0..]);

    polyz_unpack(a, result[0..]);
}

pub fn poly_challenge(c: *poly, seed: [params.CTILDEBYTES]u8) void {
    var ctx = symmetric.init256();

    var state = &ctx.state;
    state.update(seed);
    state.finalize();
    var buf: [SHAKE256_RATE]u8 = undefined;
    state.squeezeBlocks(buf[0..], 1);

    var signs: u64 = 0;
    var i: usize = 0;
    while (i < 8) : (i += 1) {
        signs |= @as(u64, @intCast(buf[i])) << (i * 8);
    }
    var pos: usize = 8;

    for (c.coeffs) |*coeff| {
        coeff.* = 0;
    }

    i = params.N - params.TAU;
    while (i < params.N) : (i += 1) {
        const b = blk: {
            while (true) {
                if (pos >= SHAKE256_RATE) {
                    state.squeeze(buf[0..]);
                    pos = 0;
                }

                const candidate = buf[pos];
                pos += 1;

                if (candidate <= @as(u16, @intCast(i))) {
                    break :blk candidate;
                }
            }
        };

        std.mem.swap(i16, &c.coeffs[i], &c.coeffs[b]);
        c.coeffs[b] = @as(i16, 1 - 2 * @as(i1, @intFromBool(signs & 1 == 1)));
        signs >>= 1;
    }
}

pub fn polyeta_pack(r: []u8, a: *const poly) void {
    var t_shrunk: [8]u8 = undefined;
    var i: usize = 0;
    if (params.ETA == 2) {
        while (i < params.N / 8) : (i += 1) {
            var j = 0;
            while (j < 8) : (j += 1) {
                t_shrunk[j] = 2 - a.coeffs[i * 8 + j];
            }

            r[i * 3 + 0] = (t_shrunk[0] >> 0) | (t_shrunk[1] << 3) | (t_shrunk[2] << 6);
            r[i * 3 + 1] = (t_shrunk[2] >> 2) | (t_shrunk[3] << 1) | (t_shrunk[4] << 4) | (t_shrunk[5] << 7);
            r[i * 3 + 2] = (t_shrunk[5] >> 1) | (t_shrunk[6] << 2) | (t_shrunk[7] << 5);
        }
    } else if (params.ETA == 4) {
        while (i < params.N / 2) : (i += 1) {
            const a0 = a.coeffs[i * 2 + 0];
            const a1 = a.coeffs[i * 2 + 1];

            const t0 = 4 - a0;
            const t1 = 4 - a1;

            r[i] = @as(u8, @truncate((t0) | (t1) << 4));
        }
    } else {
        @compileError("Unsupported value for ETA");
    }
}

pub fn polyeta_unpack(a: *poly, r: []const u8) void {
    var i: usize = 0;
    if (params.ETA == 2) {
        while (i < params.N / 8) : (i += 1) {
            var j = 0;
            while (j < 8) : (j += 1) {
                if (j == 2) {
                    a.coeffs[i * 8 + j] = (r[i * 3 + j] >> 6) | (a[3 * i + 1] << 2) & 7;
                } else if (j == 5) {
                    a.coeffs[i * 8 + j] = (r[i * 3 + j] >> 7) | (a[3 * i + 2] << 1) & 7;
                } else {
                    a.coeffs[i * 8 + j] = (r[i * 3 + j] >> (j % 8)) & 7;
                }
            }
        }
    } else if (params.ETA == 4) {
        while (i < params.N / 2) : (i += 1) {
            a.coeffs[2 * i + 0] = r[i] & 0x0F;
            a.coeffs[2 * i + 1] = r[i] >> 4;
            a.coeffs[2 * i + 0] = params.ETA - a.coeffs[2 * i + 0];
            a.coeffs[2 * i + 1] = params.ETA - a.coeffs[2 * i + 1];
        }
    } else {
        @compileError("Unsupported value for ETA");
    }
}

pub fn polyt1_pack(r: []u8, a: *const poly) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        r[5 * i + 0] = a.coeffs[4 * i + 0] >> 0;
        r[5 * i + 1] = (a.coeffs[4 * i + 0] >> 8) | (a.coeffs[4 * i + 1] << 2);
        r[5 * i + 2] = (a.coeffs[4 * i + 1] >> 6) | (a.coeffs[4 * i + 2] << 4);
        r[5 * i + 3] = (a.coeffs[4 * i + 1] >> 4) | (a.coeffs[4 * i + 2] << 6);
        r[5 * i + 4] = a.coeffs[4 * i + 2] >> 2;
    }
}

pub fn polyt1_unpack(a: *poly, r: []const u8) void {
    var i: usize = 0;
    while (i < params.N) : (i += 1) {
        a.coeffs[4 * i + 0] = r[5 * i + 0] >> 0 | (@as(u32, r[5 * i + 1]) << 8) & 0x3FF;
        a.coeffs[4 * i + 1] = r[5 * i + 1] >> 2 | (@as(u32, r[5 * i + 2]) << 6) & 0x3FF;
        a.coeffs[4 * i + 2] = r[5 * i + 2] >> 4 | (@as(u32, r[5 * i + 3]) << 4) & 0x3FF;
        a.coeffs[4 * i + 3] = r[5 * i + 3] >> 6 | (@as(u32, r[5 * i + 4]) << 2) & 0x3FF;
    }
}

pub fn polyt0_pack(r: []u8, a: *const poly) void {
    var i = 0;
    var t: [8]u32 = undefined;

    while (i < params.N / 8) : (i += 1) {
        t[0] = (1 << (params.D - 1)) - a.coeffs[8 * i + 0];
        t[1] = (1 << (params.D - 1)) - a.coeffs[8 * i + 1];
        t[2] = (1 << (params.D - 1)) - a.coeffs[8 * i + 2];
        t[3] = (1 << (params.D - 1)) - a.coeffs[8 * i + 3];
        t[4] = (1 << (params.D - 1)) - a.coeffs[8 * i + 4];
        t[5] = (1 << (params.D - 1)) - a.coeffs[8 * i + 5];
        t[6] = (1 << (params.D - 1)) - a.coeffs[8 * i + 6];
        t[7] = (1 << (params.D - 1)) - a.coeffs[8 * i + 7];

        r[13 * i + 0] = t[0];
        r[13 * i + 1] = t[0] >> 8;
        r[13 * i + 1] |= t[1] << 5;
        r[13 * i + 2] = t[1] >> 3;
        r[13 * i + 3] = t[1] >> 11;
        r[13 * i + 3] |= t[2] << 2;
        r[13 * i + 4] = t[2] >> 6;
        r[13 * i + 4] |= t[3] << 7;
        r[13 * i + 5] = t[3] >> 1;
        r[13 * i + 6] = t[3] >> 9;
        r[13 * i + 6] |= t[4] << 4;
        r[13 * i + 7] = t[4] >> 4;
        r[13 * i + 8] = t[4] >> 12;
        r[13 * i + 8] |= t[5] << 1;
        r[13 * i + 9] = t[5] >> 7;
        r[13 * i + 9] |= t[6] << 6;
        r[13 * i + 10] = t[6] >> 2;
        r[13 * i + 11] = t[6] >> 10;
        r[13 * i + 11] |= t[7] << 3;
        r[13 * i + 12] = t[7] >> 5;
    }
}

pub fn polyt0_unpack(r: *poly, a: []const u8) void {
    var i = 0;
    while (i < params.N / 8) : (i += 1) {
        r.coeffs[8 * i + 0] = a[13 * i + 0];
        r.coeffs[8 * i + 0] |= @as(u32, a[13 * i + 1]) << 8;
        r.coeffs[8 * i + 0] &= 0x1FFF;

        r.coeffs[8 * i + 1] = a[13 * i + 1] >> 5;
        r.coeffs[8 * i + 1] |= @as(u32, a[13 * i + 2]) << 3;
        r.coeffs[8 * i + 1] |= @as(u32, a[13 * i + 3]) << 11;
        r.coeffs[8 * i + 1] &= 0x1FFF;

        r.coeffs[8 * i + 2] = a[13 * i + 3] >> 2;
        r.coeffs[8 * i + 2] |= @as(u32, a[13 * i + 4]) << 6;
        r.coeffs[8 * i + 2] &= 0x1FFF;

        r.coeffs[8 * i + 3] = a[13 * i + 4] >> 7;
        r.coeffs[8 * i + 3] |= @as(u32, a[13 * i + 5]) << 1;
        r.coeffs[8 * i + 3] |= @as(u32, a[13 * i + 6]) << 9;
        r.coeffs[8 * i + 3] &= 0x1FFF;

        r.coeffs[8 * i + 4] = a[13 * i + 6] >> 4;
        r.coeffs[8 * i + 4] |= @as(u32, a[13 * i + 7]) << 4;
        r.coeffs[8 * i + 4] |= @as(u32, a[13 * i + 8]) << 12;
        r.coeffs[8 * i + 4] &= 0x1FFF;

        r.coeffs[8 * i + 5] = a[13 * i + 8] >> 1;
        r.coeffs[8 * i + 5] |= @as(u32, a[13 * i + 9]) << 7;
        r.coeffs[8 * i + 5] &= 0x1FFF;

        r.coeffs[8 * i + 6] = a[13 * i + 9] >> 6;
        r.coeffs[8 * i + 6] |= @as(u32, a[13 * i + 10]) << 2;
        r.coeffs[8 * i + 6] |= @as(u32, a[13 * i + 11]) << 10;
        r.coeffs[8 * i + 6] &= 0x1FFF;

        r.coeffs[8 * i + 7] = a[13 * i + 11] >> 3;
        r.coeffs[8 * i + 7] |= @as(u32, a[13 * i + 12]) << 5;
        r.coeffs[8 * i + 7] &= 0x1FFF;

        r.coeffs[8 * i + 0] = (1 << (params.D - 1)) - r.coeffs[8 * i + 0];
        r.coeffs[8 * i + 1] = (1 << (params.D - 1)) - r.coeffs[8 * i + 1];
        r.coeffs[8 * i + 2] = (1 << (params.D - 1)) - r.coeffs[8 * i + 2];
        r.coeffs[8 * i + 3] = (1 << (params.D - 1)) - r.coeffs[8 * i + 3];
        r.coeffs[8 * i + 4] = (1 << (params.D - 1)) - r.coeffs[8 * i + 4];
        r.coeffs[8 * i + 5] = (1 << (params.D - 1)) - r.coeffs[8 * i + 5];
        r.coeffs[8 * i + 6] = (1 << (params.D - 1)) - r.coeffs[8 * i + 6];
        r.coeffs[8 * i + 7] = (1 << (params.D - 1)) - r.coeffs[8 * i + 7];
    }
}

pub fn polyz_pack(r: []u8, a: *const poly) void {
    var i = 0;
    var t: [4]u32 = undefined;

    if (params.GAMMA1 == (1 << 17)) {
        while (i < params.N / 4) : (i += 1) {
            t[0] = params.GAMMA1 - a.coeffs[4 * i + 0];
            t[1] = params.GAMMA1 - a.coeffs[4 * i + 1];
            t[2] = params.GAMMA1 - a.coeffs[4 * i + 2];
            t[3] = params.GAMMA1 - a.coeffs[4 * i + 3];
            r[9 * i + 0] = t[0];
            r[9 * i + 1] = t[0] >> 8;
            r[9 * i + 2] = t[0] >> 16;
            r[9 * i + 2] |= t[1] << 2;
            r[9 * i + 3] = t[1] >> 6;
            r[9 * i + 4] = t[1] >> 14;
            r[9 * i + 4] |= t[2] << 4;
            r[9 * i + 5] = t[2] >> 4;
            r[9 * i + 6] = t[2] >> 12;
            r[9 * i + 6] |= t[3] << 6;
            r[9 * i + 7] = t[3] >> 2;
            r[9 * i + 8] = t[3] >> 10;
        }
    } else if (params.GAMMA1 == (1 << 19)) {
        while (i < params.N / 2) : (i += 1) {
            t[0] = params.GAMMA1 - a.coeffs[2 * i + 0];
            t[1] = params.GAMMA1 - a.coeffs[2 * i + 1];
            r[5 * i + 0] = t[0];
            r[5 * i + 1] = t[0] >> 8;
            r[5 * i + 2] = t[0] >> 16;
            r[5 * i + 2] |= t[1] << 4;
            r[5 * i + 3] = t[1] >> 4;
            r[5 * i + 4] = t[1] >> 12;
        }
    }
}

pub fn polyz_unpack(r: *poly, a: []const u8) void {
    var i: usize = 0;
    if (params.GAMMA1 == (1 << 17)) {
        while (i < params.N / 4) : (i += 1) {
            r.coeffs[4 * i + 0] = a[9 * i + 0];
            r.coeffs[4 * i + 0] |= @as(u32, a[9 * i + 1]) << 8;
            r.coeffs[4 * i + 0] |= @as(u32, a[9 * i + 2]) << 16;
            r.coeffs[4 * i + 0] &= 0x3FFFF;

            r.coeffs[4 * i + 1] = a[9 * i + 2] >> 2;
            r.coeffs[4 * i + 1] |= @as(u32, a[9 * i + 3]) << 6;
            r.coeffs[4 * i + 1] |= @as(u32, a[9 * i + 4]) << 14;
            r.coeffs[4 * i + 1] &= 0x3FFFF;

            r.coeffs[4 * i + 2] = a[9 * i + 4] >> 4;
            r.coeffs[4 * i + 2] |= @as(u32, a[9 * i + 5]) << 4;
            r.coeffs[4 * i + 2] |= @as(u32, a[9 * i + 6]) << 12;
            r.coeffs[4 * i + 2] &= 0x3FFFF;

            r.coeffs[4 * i + 3] = a[9 * i + 6] >> 6;
            r.coeffs[4 * i + 3] |= @as(u32, a[9 * i + 7]) << 2;
            r.coeffs[4 * i + 3] |= @as(u32, a[9 * i + 8]) << 10;
            r.coeffs[4 * i + 3] &= 0x3FFFF;

            r.coeffs[4 * i + 0] = params.GAMMA1 - r.coeffs[4 * i + 0];
            r.coeffs[4 * i + 1] = params.GAMMA1 - r.coeffs[4 * i + 1];
            r.coeffs[4 * i + 2] = params.GAMMA1 - r.coeffs[4 * i + 2];
            r.coeffs[4 * i + 3] = params.GAMMA1 - r.coeffs[4 * i + 3];
        }
    } else if (params.GAMMA1 == (1 << 19)) {
        while (i < params.N / 2) : (i += 1) {
            r.coeffs[2 * i + 0] = a[5 * i + 0];
            r.coeffs[2 * i + 0] |= @as(u32, a[5 * i + 1]) << 8;
            r.coeffs[2 * i + 0] |= @as(u32, a[5 * i + 2]) << 16;
            r.coeffs[2 * i + 0] &= 0xFFFFF;

            r.coeffs[2 * i + 1] = a[5 * i + 2] >> 4;
            r.coeffs[2 * i + 1] |= @as(u32, a[5 * i + 3]) << 4;
            r.coeffs[2 * i + 1] |= @as(u32, a[5 * i + 4]) << 12;
            r.coeffs[2 * i + 1] &= 0xFFFFF;
            r.coeffs[2 * i + 0] = params.GAMMA1 - r.coeffs[2 * i + 0];
            r.coeffs[2 * i + 1] = params.GAMMA1 - r.coeffs[2 * i + 1];
        }
    }
}
