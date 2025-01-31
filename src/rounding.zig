const params = @import("params.zig");

const D = params.D;
const Q = params.Q;
const GAMMA2 = params.GAMMA2;

pub fn power2round(a0: *i32, a: i32) i32 {
    const a1: i32 = (a + (1 << (D - 1)) - 1) >> D;
    a0 = a - (a1 << D);
    return a1;
}

pub fn decompose(a0: *i32, a: i32) i32 {
    var a1 = (a + 127) >> 7;
    switch (GAMMA2) {
        (Q - 1) / 32 => {
            a1 = (a1 * 1025 + (1 << 21)) >> 22;
            a1 &= 15;
        },
        (Q - 1) / 88 => {
            a1 = (a1 * 11275 + (1 << 23)) >> 24;
            a1 ^= ((43 - a1) >> 31) & a1;
        },
        _ => unreachable,
    }

    a0 = a - a1 * 2 * GAMMA2;
    a0 -= (((Q - 1) / 2 - *a0) >> 31) & Q;
    return a1;
}

pub fn make_hint(a0: i32, a1: i32) u8 {
    switch (a0) {
        a0 > GAMMA2 => return 1,
        a0 < -GAMMA2 => return 1,
        a0 == -GAMMA2 => if (a1 != 0) return 1 else return 0,
        else => return 0,
    }
}

pub fn use_hint(a: i32, hint: u8) i32 {
    var a0: i32 = undefined;

    const a1 = decompose(&a0, a);
    if (hint == 0) return a0;
    switch (GAMMA2) {
        (Q - 1) / 32 => {
            if (a0 > 0) return (a1 + 1) & 15 else return (a1 - 1) & 15;
        },
        (Q - 1) / 88 => {
            if (a0 > 0) {
                return if (a1 == 43) return 0 else a1 + 1;
            } else {
                return if (a1 == 0) return 43 else a1 - 1;
            }
        },
        _ => unreachable,
    }
}
