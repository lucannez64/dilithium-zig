const std = @import("std");
const params = @import("params.zig");

const Shake128 = std.crypto.hash.sha3.Shake128;
const Shake256 = std.crypto.hash.sha3.Shake256;

pub const Shake128_state = struct {
    state: Shake128,
    s: [25]u64,
};

pub const Shake256_state = struct {
    state: Shake256,
    s: [25]u64,
};

pub fn init128() Shake128_state {
    var state: Shake128_state = undefined;
    state.state = Shake128.init(.{});
    state.s = [_]u64{0} ** 25;
    return state;
}

pub fn init256() Shake256_state {
    var state: Shake256_state = undefined;
    state.s = [_]u64{0} ** 25;
    state.state = Shake256.init(.{});
    return state;
}

pub fn dilithium_shake128_stream_init(state: *Shake128_state, seed: [params.SEEDBYTES]u8, nonce: u16) void {
    var t: [2]u8 = undefined;
    t[0] = @truncate(nonce);
    t[1] = @truncate(nonce >> 8);
    state = Shake128.init(.{}) catch unreachable;
    Shake128.update(state.state, &seed);
    Shake128.update(state.state, &t);
    Shake128.final(state.state, state.s);
}

pub fn dilithium_shake256_stream_init(state: *Shake256_state, seed: [params.CRHBYTES]u8, nonce: u16) void {
    var t: [2]u8 = undefined;
    t[0] = @truncate(nonce);
    t[1] = @truncate(nonce >> 8);
    state = Shake128.init(.{}) catch unreachable;
    Shake128.update(state.state, &seed);
    Shake128.update(state.state, &t);
    Shake128.final(state.state, state.s);
}

pub fn absorb(state: *Shake256_state, data: []const u8) void {
    state.state.st.absorb(data);
}

pub fn squeeze(state: *Shake256_state, data: []u8) void {
    state.state.st.squeeze(data);
}
