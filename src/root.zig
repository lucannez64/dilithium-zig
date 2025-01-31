const std = @import("std");
pub const params = @import("params.zig");
const sign = @import("sign.zig");
const testing = std.testing;

pub fn dilithium5_keypair(pk: []u8, sk: []u8) void {
    sign.crypto_sign_keypair(&pk, &sk);
}

pub fn dilithium5_sign(sk: []const u8, msg: []const u8, sig: []u8, ctx: []const u8) sign.SignatureError!void {
    try sign.crypto_sign_signature(&sig, msg, ctx, sk);
}

pub fn dilithium5_verify(pk: []const u8, msg: []const u8, sig: []const u8, ctx: []const u8) sign.SignatureError!bool {
    return sign.crypto_sign_verify(&sig, msg, ctx, pk);
}

pub fn dilithium5_open(msg: []u8, sig: []const u8, ctx: []const u8, pk: []const u8) sign.SignatureError!void {
    return sign.crypto_sign_open(msg, sig, ctx, pk);
}

test "basic add functionality" {}
