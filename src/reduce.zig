const params = @import("params.zig");
const Q = params.Q;
const std = @import("std");

const MONT = -4186625; // 2^32 % Q
const QINV = 58728449; // q^(-1) mod 2^32

pub fn montgomery_reduce(a: i64) i32 {
    const t: i32 = @truncate(a * QINV);
    return @truncate((a - @as(i64, t) * Q) >> 32);
}

pub fn reduce32(a: i32) i32 {
    return a - ((a + (1 << 22)) >> 23) * Q;
}

pub fn caddq(a: i32) i32 {
    return a + ((a >> 31) & Q);
}

pub fn freeze(a: i32) i32 {
    return caddq(reduce32(a));
}
